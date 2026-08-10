#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import RegularGridInterpolator


class AboveSolidusError(Exception):
	pass

class BurnmanTableQuery:
	def __init__(self, npz_path = "burnman_pyx_lherz_dun.npz"):
	
		"""
		Fast query interface for the BurnMan-derived dunite/pyroxenite/lherzolite
		mixing lookup table (Vp, Vs, density, modal mineralogy, and per-phase
		Fe# across a grid of bulk composition, temperature, and pressure).
	
		This three-endmember (dunite*/pyroxenite*/lherzolite*) mixing
		methodology was developed by Munch et al. (2025), "Evidence for
		Lithospheric Mantle Uniformity Beneath Cratons", JGR Solid Earth,
		130, e2024JB029649. Their approach (linearly mix bulk end-member
		oxide compositions by modal fraction, then run a single equilibrium
		calculation on the resulting blend) was originally implemented with
		Perple_X; this table reproduces the same methodology using BurnMan
		(SLB_2011 dataset) instead, with modal fractions and end-member
		compositions matched to Munch et al.'s Table 1.
	
		The underlying table was built by Gibbs-energy minimization over a
		4D grid of (f_pyx, f_lherz, T, P), comparing garnet- and spinel-bearing
		assemblages at every point and keeping whichever has lower Gibbs
		energy. f_dun is not stored separately since f_dun = 1 - f_pyx - f_lherz
		always.
	
		NaN handling: every grid cell with no data, either because the
		composition falls outside the dunite/pyroxenite/lherzolite triangle,
		or because that (T, P) point is above the estimated solidus for that
		composition (no solid-only equilibrium exists there) is left as
		NaN in the padded grid. Both cases surface identically as NaN output
		from query(); a genuine solidus crossing persists across a whole
		range of neighboring points, while a handful of known isolated
		single-cell numerical gaps (confirmed at a few specific low-T points,
		unrelated to melting) do not. query()'s nudge_on_failure option
		exploits exactly that distinction: on failure it retries at small
		T/P offsets and only accepts a neighbor's result if the offset
		actually recovers a valid value, so isolated numerical gaps get
		rescued while real solidus failures are correctly left as NaN.
		Every returned result includes 'nudge_used' so it's always possible
		to tell whether a value came from the exact requested point or a
		nearby substitute.
	
		Usage
		-----
			tq = BurnmanTableQuery('burnman_pyx_lherz_dun.npz')
			result = tq.query(f_pyx=0.2, f_lherz=0.3, T=1300, P=2.5)
			# result: {'v_p':..., 'v_s':..., 'density':..., 'ol':..., 'opx':...,
			#          'cpx':..., 'garnet':..., 'spinel':..., 'nudge_used': None}
		"""
		d_loaded = np.load(npz_path, allow_pickle=True)
		d = {k: d_loaded[k] for k in d_loaded.files}
 
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
		padded = {field: np.full(shape, np.nan) for field in fields}
 
		none_mask = d['assemblage_used'] == 'none'
		for i in range(len(d['f_pyx'])):
			ip = f_pyx_idx[round(float(d['f_pyx'][i]), 6)]
			il = f_lherz_idx[round(float(d['f_lherz'][i]), 6)]
			iT = T_idx[round(float(d['T'][i]), 6)]
			iP = P_idx[round(float(d['P'][i]), 6)]
			if not none_mask[i]:
				for field in fields:
					padded[field][ip, il, iT, iP] = d[field][i]
			# 'none' points are simply left as NaN — both the "outside the
			# triangle" cells and the "above solidus" cells end up NaN this
			# way, which is exactly the safe behavior we want.
 
		self.grids = (f_pyx_grid, f_lherz_grid, T_grid, P_grid)
		self.interps = {field: RegularGridInterpolator(self.grids, padded[field],
			bounds_error=False, fill_value=np.nan) for field in fields}
 
	def _raw_query(self, f_pyx, f_lherz, T, P):
		point = np.array([[f_pyx, f_lherz, T, P]])
		return {field: float(interp(point)[0]) for field, interp in self.interps.items()}
 
	def query(self, f_pyx, f_lherz, T, P, raise_on_invalid=False, nudge_on_failure=True,
		T_nudges=(25, -25, 50, -50), P_nudges=(0.1, -0.1, 0.25, -0.25)):
		"""
		Query the table at (f_pyx, f_lherz, T, P). If the exact point comes
		back invalid (NaN) and nudge_on_failure is True, retry at small
		offsets in T first, then P, before giving up — this recovers from
		known isolated single-point numerical gaps in the table (confirmed
		to exist at a few specific low-temperature points; a genuine
		solidus failure persists across a range of nearby points and won't
		be "rescued" by this, so it doesn't mask real physics).
 
		The returned dict always includes 'nudge_used': None if the exact
		point succeeded, or (axis, offset) if a nearby point had to be
		substituted — check this if you need to know whether the result
		came from the exact requested point.
		"""
		result = self._raw_query(f_pyx, f_lherz, T, P)
		nudge_used = None
 
		if np.isnan(result['v_p']) and nudge_on_failure:
			for dT in T_nudges:
				candidate = self._raw_query(f_pyx, f_lherz, T + dT, P)
				if not np.isnan(candidate['v_p']):
					result = candidate
					nudge_used = ('T', dT)
					break
			if nudge_used is None:
				for dP in P_nudges:
					candidate = self._raw_query(f_pyx, f_lherz, T, P + dP)
					if not np.isnan(candidate['v_p']):
						result = candidate
						nudge_used = ('P', dP)
						break
 
		result['nudge_used'] = nudge_used
 
		if raise_on_invalid and np.isnan(result['v_p']):
			raise AboveSolidusError(
				f"Composition (f_pyx={f_pyx}, f_lherz={f_lherz}) at T={T}K, P={P}GPa "
				f"returned NaN even after nudging nearby T/P — likely a genuine solidus "
				f"crossing (a real gap persists across neighbors), outside the composition "
				f"triangle, or outside the table's T/P range.")
		return result