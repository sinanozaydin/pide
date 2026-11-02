#!/usr/bin/env python3

import os
import warnings

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import distance_transform_edt


class TableLookupError(Exception):
	pass


#older name, kept so existing imports and except clauses keep working
AboveSolidusError = TableLookupError


class TriangleTableQuery:
	def __init__(self, npz_path = None, warn_missing = True):

		"""
		Fast query interface for the dunite/pyroxenite/lherzolite mixing
		lookup table (Vp, Vs, density, modal mineralogy, and per-phase Fe#
		across a grid of bulk composition, temperature, and pressure).

		This three-endmember (dunite*/pyroxenite*/lherzolite*) mixing
		methodology was developed by Munch et al. (2025), "Evidence for
		Lithospheric Mantle Uniformity Beneath Cratons", JGR Solid Earth,
		130, e2024JB029649: linearly mix bulk end-member oxide compositions
		by modal fraction, then run a single equilibrium calculation on the
		resulting blend. End-member compositions and modal fractions here
		are matched to their Table 1.

		This class is agnostic about which code produced the table. It
		needs a flat npz carrying f_pyx, f_lherz, T, P, assemblage_used and
		the fields listed below, one row per grid point, on a regular
		4D grid of (f_pyx, f_lherz, T, P). f_dun is not stored separately
		since f_dun = 1 - f_pyx - f_lherz always.

		NaN handling. A cell is NaN where the table records no assemblage,
		which means the equilibrium calculation returned nothing there.
		Compositions outside the triangle are rejected before the
		interpolator sees them, so they are NaN for a stated reason rather
		than by accident.

		nudge_on_failure exploits the geometry of the remaining gaps. A gap
		with a physical cause persists across a range of neighbouring
		points; an isolated numerical failure does not. On failure, query()
		retries at small T then P offsets and accepts a neighbour's result
		only if the offset actually recovers a value, so single-cell gaps
		are rescued while broad regions are correctly left as NaN. Every
		result carries 'nudge_used' so you can always tell whether a value
		came from the requested point or a substitute.

		Usage
		-----
			tq = TriangleTableQuery('magemin_pyx_lherz_dun.npz')
			result = tq.query(f_pyx=0.2, f_lherz=0.3, T=1300, P=2.5)
			# result: {'v_p':..., 'v_s':..., 'density':..., 'ol':..., 'opx':...,
			#          'cpx':..., 'garnet':..., 'spinel':..., 'nudge_used': None}
		"""

		if npz_path is None:
			here = os.path.dirname(os.path.abspath(__file__))
			#the current default first, then older names, so an existing
			#installation keeps working without being renamed
			for name in ('pyx_lherz_dun_table.npz',
				'magemin_pyx_lherz_dun.npz'):
				candidate = os.path.join(here, name)
				if os.path.exists(candidate):
					npz_path = candidate
					break
			else:
				raise FileNotFoundError(
					'no lookup table beside %s. Looked for '
					'pyx_lherz_dun_table.npz and '
					'magemin_pyx_lherz_dun.npz. Pass npz_path '
					'explicitly.' % here)

		d_loaded = np.load(npz_path, allow_pickle=True)
		d = {k: d_loaded[k] for k in d_loaded.files}
		self.npz_path = npz_path

		f_pyx_grid = np.array(sorted(set(np.round(d['f_pyx'], 6))))
		f_lherz_grid = np.array(sorted(set(np.round(d['f_lherz'], 6))))
		T_grid = np.array(sorted(set(np.round(d['T'], 6))))
		P_grid = np.array(sorted(set(np.round(d['P'], 6))))

		f_pyx_idx = {v: i for i, v in enumerate(f_pyx_grid)}
		f_lherz_idx = {v: i for i, v in enumerate(f_lherz_grid)}
		T_idx = {v: i for i, v in enumerate(T_grid)}
		P_idx = {v: i for i, v in enumerate(P_grid)}

		shape = (len(f_pyx_grid), len(f_lherz_grid), len(T_grid), len(P_grid))

		fields = ['v_p', 'v_s', 'density', 'ol', 'opx', 'cpx', 'garnet', 'spinel',
			'ol_xfe', 'opx_xfe', 'cpx_xfe', 'garnet_xfe', 'spinel_xfe',
			'SiO2', 'MgO', 'FeO', 'Al2O3', 'CaO', 'Na2O']

		#a field the npz lacks would otherwise become an all-NaN
		#interpolator that returns NaN for every query forever, with
		#nothing to say why
		missing = [f for f in fields if f not in d]
		if missing and warn_missing:
			warnings.warn('%s has no %s; those fields will be NaN for '
				'every query' % (os.path.basename(npz_path),
				', '.join(missing)), RuntimeWarning, stacklevel=2)
		fields = [f for f in fields if f in d]
		self.fields = list(fields)
		self.missing_fields = missing

		padded = {field: np.full(shape, np.nan) for field in fields}

		none_mask = np.asarray(d['assemblage_used']) == 'none'
		inside = np.zeros(shape[:2], bool)
		for i in range(len(d['f_pyx'])):
			ip = f_pyx_idx[round(float(d['f_pyx'][i]), 6)]
			il = f_lherz_idx[round(float(d['f_lherz'][i]), 6)]
			iT = T_idx[round(float(d['T'][i]), 6)]
			iP = P_idx[round(float(d['P'][i]), 6)]
			#the composition exists in the table whether or not this
			#particular (T, P) converged
			inside[ip, il] = True
			if not none_mask[i]:
				for field in fields:
					padded[field][ip, il, iT, iP] = d[field][i]
			# 'none' points are left as NaN. Those are genuine holes in
			# the equilibrium calculation and must stay visible.

		# The triangle is stored in a dense rectangle whose corner lies
		# outside it. RegularGridInterpolator's linear stencil spans 2x2
		# in (f_pyx, f_lherz), so any query whose cell touches an
		# out-of-triangle node returns NaN, because 0 * NaN is NaN even
		# where the weight is zero. Left alone that makes every
		# composition with f_dun below one grid step unreachable: 41 of
		# the 231 nodes at f_step 0.05, an invisible wall a sampler keeps
		# hitting and counting as a rejection.
		#
		# Fill that corner by nearest neighbour from inside the triangle.
		# Inside cells index themselves, so real values and genuine
		# 'none' holes are both untouched, and the filler only ever
		# carries weight for queries sitting on the boundary. Queries
		# genuinely outside the triangle are rejected in _raw_query, so
		# the filler is never reachable as an answer in its own right.
		self._n_filled = int((~inside).sum())
		if self._n_filled:
			_, (si, sj) = distance_transform_edt(~inside,
				return_indices=True)
			for field in fields:
				padded[field] = padded[field][si, sj]

		self.grids = (f_pyx_grid, f_lherz_grid, T_grid, P_grid)
		self.interps = {field: RegularGridInterpolator(self.grids, padded[field],
			bounds_error=False, fill_value=np.nan) for field in fields}

	def in_triangle(self, f_pyx, f_lherz, tol = 1e-9):
		"""Whether this composition exists at all. f_dun = 1 - f_pyx -
		f_lherz cannot be negative and neither fraction can be."""
		return (f_pyx >= -tol and f_lherz >= -tol
			and f_pyx + f_lherz <= 1.0 + tol)

	def _nan_result(self):
		return {field: np.nan for field in self.interps}

	def _raw_query(self, f_pyx, f_lherz, T, P):
		#outside the triangle there is no rock, and the padding fill
		#would otherwise answer with a neighbour's value
		if not self.in_triangle(f_pyx, f_lherz):
			return self._nan_result()
		point = np.array([[f_pyx, f_lherz, T, P]])
		return {field: float(interp(point)[0])
			for field, interp in self.interps.items()}

	@staticmethod
	def _invalid(result):
		"""A result is invalid if ANY velocity is NaN, not only v_p."""
		for k in ('v_p', 'v_s'):
			v = result.get(k, np.nan)
			if v is None or np.isnan(v):
				return True
		return False

	def query(self, f_pyx, f_lherz, T, P, raise_on_invalid=False, nudge_on_failure=True,
		T_nudges=(25, -25, 50, -50), P_nudges=(0.1, -0.1, 0.25, -0.25)):
		"""
		Query the table at (f_pyx, f_lherz, T, P). If the exact point comes
		back invalid (NaN) and nudge_on_failure is True, retry at small
		offsets in T first, then P, before giving up. This recovers from
		isolated single-point gaps in the table. A gap with a physical
		cause persists across a range of nearby points and will not be
		"rescued" by this, so it does not mask real behaviour.

		The returned dict always includes 'nudge_used': None if the exact
		point succeeded, or (axis, offset) if a nearby point had to be
		substituted. Check it if you need to know whether the result came
		from the exact requested point.
		"""
		outside = not self.in_triangle(f_pyx, f_lherz)
		result = self._raw_query(f_pyx, f_lherz, T, P)
		nudge_used = None

		#nudging in T and P cannot rescue a composition that does not
		#exist, so do not spend eight interpolations discovering that
		if self._invalid(result) and nudge_on_failure and not outside:
			for dT in T_nudges:
				candidate = self._raw_query(f_pyx, f_lherz, T + dT, P)
				if not self._invalid(candidate):
					result = candidate
					nudge_used = ('T', dT)
					break
			if nudge_used is None:
				for dP in P_nudges:
					candidate = self._raw_query(f_pyx, f_lherz, T, P + dP)
					if not self._invalid(candidate):
						result = candidate
						nudge_used = ('P', dP)
						break

		result['nudge_used'] = nudge_used

		if raise_on_invalid and self._invalid(result):
			if outside:
				raise TableLookupError(
					f"Composition (f_pyx={f_pyx}, f_lherz={f_lherz}) is "
					f"outside the triangle: f_dun = "
					f"{1.0 - f_pyx - f_lherz:.4f} cannot be negative.")
			raise TableLookupError(
				f"Composition (f_pyx={f_pyx}, f_lherz={f_lherz}) at T={T}K, P={P}GPa "
				f"returned NaN even after nudging nearby T/P. The gap persists "
				f"across neighbouring points, so it is either a region the table "
				f"has no solution for, or outside the table's T/P range.")
		return result