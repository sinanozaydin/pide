import numpy as np
import burnman
from burnman import equilibrate
from burnman.minerals import SLB_2011

"""
Registry of named reference mantle bulk compositions, for building separate
KD lookup tables (PUM, DMM, etc.) that pide users can choose between.

Add new compositions here as {'name': (SiO2, MgO, FeO, Al2O3, CaO, Na2O)}
in wt% oxides. Conversion to cation moles (with exact charge balance) is
handled by to_cations() below, do not hand-round the converted values,
BurnMan enforces charge balance to a tolerance of 1e-12.
"""

_OXIDE_MOLAR_MASS = {
	'SiO2': 60.08, 'MgO': 40.30, 'FeO': 71.85,
	'Al2O3': 101.96, 'CaO': 56.08, 'Na2O': 61.98,
}

# wt% oxides: (SiO2, MgO, FeO, Al2O3, CaO, Na2O)
COMPOSITIONS_WT_PCT = {
	'PUM': (45.00, 37.80, 8.05, 4.45, 3.55, 0.36),   # McDonough & Sun (1995), Primitive Upper Mantle
	'DMM': (44.71, 38.73, 8.18, 3.98, 3.17, 0.13),   # Workman & Hart (2005), Depleted MORB Mantle
	# Add further reference compositions here, e.g. 'KLB-1': (...),
	# only add values you've verified against the source paper — see the
	# Mg# cross-check pattern used for PUM and DMM above before trusting a
	# new entry.
}

def to_cations(SiO2, MgO, FeO, Al2O3, CaO, Na2O):
	"""Convert wt% oxides to cation mole amounts + charge-balanced O,
	in the format burnman.equilibrate expects."""
	Si = SiO2 / _OXIDE_MOLAR_MASS['SiO2']
	Mg = MgO / _OXIDE_MOLAR_MASS['MgO']
	Fe = FeO / _OXIDE_MOLAR_MASS['FeO']
	Al = 2 * Al2O3 / _OXIDE_MOLAR_MASS['Al2O3']
	Ca = CaO / _OXIDE_MOLAR_MASS['CaO']
	Na = 2 * Na2O / _OXIDE_MOLAR_MASS['Na2O']
	O = 0.5 * Na + Ca + Mg + Fe + 1.5 * Al + 2 * Si
	return {'Na': Na, 'Fe': Fe, 'Mg': Mg, 'Si': Si, 'Ca': Ca, 'Al': Al, 'O': O}


def get_composition(name):
	"""Return the cation-mole composition dict for a named reference composition."""
	if name not in COMPOSITIONS_WT_PCT:
		raise ValueError(f"Unknown composition '{name}'. Available: {list(COMPOSITIONS_WT_PCT)}")
	return to_cations(*COMPOSITIONS_WT_PCT[name])
	
_DEFAULT_T_GRID = np.linspace(800, 2000, 25)
_DEFAULT_P_GRID = np.linspace(1, 8.0, 25)
_DEFAULT_XFE_GRID = np.linspace(0.03, 0.25, 20)

def _new_assemblage():
	"""Fresh olivine+opx+cpx+garnet phase objects and starting guesses.
	Called per (T,P) column so each column starts from the same generic
	guess rather than carrying over a possibly-bad state."""
	
	ol = SLB_2011.mg_fe_olivine()
	opx = SLB_2011.orthopyroxene()
	cpx = SLB_2011.clinopyroxene()
	gt = SLB_2011.garnet()
	
	assemblage = burnman.Composite(phases=[ol, opx, cpx, gt], fractions=[0.6, 0.15, 0.15, 0.1])
	
	ol.set_composition([0.90, 0.10])
	opx.set_composition([0.80, 0.08, 0.05, 0.07])
	cpx.set_composition([0.75, 0.05, 0.02, 0.08, 0.10])
	gt.set_composition([0.65, 0.10, 0.20, 0.02, 0.03])
	
	return ol, opx, cpx, gt, assemblage

def _femg(x):
	return x / (1 - x)

def build_KD_table(composition_name=None, custom_composition_wt_pct=None, custom_name='custom',
	T_grid=None, P_grid=None, xfe_grid=None, output_path=None):
	"""
	Build a single K_D(T, P, bulk_xfe) lookup table, either for a named
	reference composition already registered in mantle_compositions.py,
	OR for a manually entered custom bulk composition.

	Parameters
	----------
	composition_name : str, optional
		Name registered in COMPOSITIONS_WT_PCT or COMPOSITIONS_MOLAR_OXIDE
		(e.g. 'PUM', 'DMM', 'KLB-1'). Mutually exclusive with
		custom_composition_wt_pct.
	custom_composition_wt_pct : dict, optional
		Manually entered bulk composition in wt% oxides, with keys
		'SiO2', 'MgO', 'FeO', 'Al2O3', 'CaO', 'Na2O'. Use this to build a
		table for a composition not in the registry (e.g. your own
		xenolith-derived or region-specific composition), without needing
		to add it to mantle_compositions.py first.
		Example: {'SiO2': 44.9, 'MgO': 38.0, 'FeO': 8.3, 'Al2O3': 4.0,
		          'CaO': 3.4, 'Na2O': 0.3}
	custom_name : str
		Label used for progress messages and the default output filename
		when using custom_composition_wt_pct. Ignored if composition_name
		is given.
	T_grid, P_grid, xfe_grid : array, optional
		As before.
	output_path : str, optional
		Defaults to 'KD_lookup_table_<composition_name or custom_name>.npz'.
	verbose : bool
		Print progress and failure counts.

	Returns
	-------
	dict with keys T_grid, P_grid, xfe_grid, KD_opx_ol_table, KD_cpx_ol_table,
	KD_gt_ol_table (also saved to output_path).
	"""
	
	if composition_name is not None and custom_composition_wt_pct is not None:
		raise ValueError("Provide either composition_name OR custom_composition_wt_pct, not both.")
	if composition_name is None and custom_composition_wt_pct is None:
		raise ValueError("Provide either composition_name or custom_composition_wt_pct.")

	if composition_name is not None:
		base_comp = get_composition(composition_name)
		label = composition_name
	else:
		required_keys = {'SiO2', 'MgO', 'FeO', 'Al2O3', 'CaO', 'Na2O'}
		missing = required_keys - set(custom_composition_wt_pct)
		if missing:
			raise ValueError(f"custom_composition_wt_pct is missing keys: {missing}")
		base_comp = to_cations(**{k: custom_composition_wt_pct[k] for k in required_keys})
		label = custom_name

	T_grid = _DEFAULT_T_GRID if T_grid is None else np.asarray(T_grid)
	P_grid = _DEFAULT_P_GRID if P_grid is None else np.asarray(P_grid)
	xfe_grid = _DEFAULT_XFE_GRID if xfe_grid is None else np.asarray(xfe_grid)
	output_path = output_path or f'KD_lookup_table_{label}.npz'

	FeMg_total = base_comp['Fe'] + base_comp['Mg']
	Si, Al, Ca, Na = base_comp['Si'], base_comp['Al'], base_comp['Ca'], base_comp['Na']

	KD_opx_ol_table = np.full((len(T_grid), len(P_grid), len(xfe_grid)), np.nan)
	KD_cpx_ol_table = np.full((len(T_grid), len(P_grid), len(xfe_grid)), np.nan)
	KD_gt_ol_table = np.full((len(T_grid), len(P_grid), len(xfe_grid)), np.nan)

	n_fail = 0
	for i, T_K in enumerate(T_grid):
		for j, P_GPa in enumerate(P_grid):
			ol, opx, cpx, gt, assemblage = _new_assemblage()

			for k, xfe in enumerate(xfe_grid):
				Fe = xfe * FeMg_total
				Mg = (1 - xfe) * FeMg_total
				O = 0.5 * Na + Ca + Mg + Fe + 1.5 * Al + 2 * Si
				composition = {'Na': Na, 'Fe': Fe, 'Mg': Mg, 'Si': Si, 'Ca': Ca, 'Al': Al, 'O': O}
				eq_constraints = [('P', P_GPa * 1e9), ('T', float(T_K))]

				try:
					sol, prm = equilibrate(composition, assemblage, eq_constraints)
				except Exception:
					sol = None

				if sol is not None and sol.success:
					ol_xfe = ol.molar_fractions[1]
					KD_opx_ol_table[i, j, k] = _femg(opx.molar_fractions[1]) / _femg(ol_xfe)
					KD_cpx_ol_table[i, j, k] = _femg(cpx.molar_fractions[1]) / _femg(ol_xfe)
					KD_gt_ol_table[i, j, k] = _femg(gt.molar_fractions[1]) / _femg(ol_xfe)
				else:
					n_fail += 1

	total_points = len(T_grid) * len(P_grid) * len(xfe_grid)

	result = dict(T_grid=T_grid, P_grid=P_grid, xfe_grid=xfe_grid,
		KD_opx_ol_table=KD_opx_ol_table, KD_cpx_ol_table=KD_cpx_ol_table,
		KD_gt_ol_table=KD_gt_ol_table)

	np.savez(output_path, **result)
	print(f"  saved to {output_path}")

	return result

