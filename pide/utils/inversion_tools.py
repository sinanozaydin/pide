import copy
import numpy as np

#Arrays larger than this are treated as static reference data and shared
#rather than copied. snapshot_report tells you whether anything mutable is
#being skipped because of it.
SNAPSHOT_MAX_ELEMENTS = 200000

SNAPSHOT_SKIP = set([
	#conductivity model database
	'cond_data_array', 'materials_data', 'mechanism_model',
	'amp_cond_data', 'amphibolite_cond_data', 'basalt_cond_data',
	'cpx_cond_data', 'fluid_cond_data', 'gabbro_cond_data',
	'garnet_cond_data', 'gneiss_cond_data', 'granite_cond_data',
	'granulite_cond_data', 'graphite_cond_data', 'kfelds_cond_data',
	'melt_cond_data', 'mica_cond_data', 'mixture_cond_data',
	'mud_cond_data', 'ol_cond_data', 'opx_cond_data',
	'other_cond_data', 'other_rock_cond_data', 'perovskite_cond_data',
	'plag_cond_data', 'quartz_cond_data', 'rwd_wds_cond_data',
	'sandstone_cond_data', 'spinel_cond_data', 'sulphides_cond_data',
	#conductivity model parameter columns
	'alpha_p', 'alpha_p_err', 'h_i', 'h_i_err',
	'h_p', 'h_p_err', 'h_pol', 'h_pol_err',
	'mg_cond', 'p_max', 'p_min', 'r', 'r_err',
	'sigma_i', 'sigma_i_err', 'sigma_p', 'sigma_p_err',
	'sigma_pol', 'sigma_pol_err', 't_max', 't_min',
	'w_calib', 'wtype', 'name', 'type',
	#melt composition reference tables
	'average_melt_composition_data', 'average_melt_composition_names',
	'melt_composition_data',
	#water partitioning and solubility tables
	'mineral_sol_calib', 'mineral_sol_fug', 'mineral_sol_name',
	'mineral_sol_o2_fug',
	'water_melt_part_function', 'water_melt_part_name',
	'water_melt_part_pchange', 'water_melt_part_type',
	'water_ol_part_function', 'water_ol_part_name',
	'water_ol_part_pchange', 'water_ol_part_type',
	'water_rwd_wds_part_function', 'water_rwd_wds_part_name',
	'water_rwd_wds_part_type',
	#reference strings and material constants
	'bib_ref', 'comp_ref', 'mat_ref', 'dens_mat'
])

#Whether lists, tuples, dicts and sets are copied.
#
#For a pide object set this to False. Every container it holds is either
#static reference data loaded from CSV at setup (the conductivity model
#tables, which dominate the snapshot cost), or a list of references to
#arrays that are snapshotted separately. mineral_frac_list and
#rock_frac_list are the second kind: they hold the very same array objects
#as self.quartz_frac and friends, and because restore_object_state writes
#arrays back IN PLACE rather than rebinding them, restoring those arrays
#restores what the lists see as well.
#
#Leave it True for any object where a container is genuinely mutated and
#is not just a view onto arrays that are already captured.
SNAPSHOT_CONTAINERS = True

#Explicit whitelist. None means "snapshot everything that is not skipped",
#which is the safe default. Set it to a set of attribute names to capture
#only those.
#
#Do not write this list from memory. Run ChangeWatcher over a few thousand
#real iterations and use the list it prints, which is what actually changed
#rather than what you remember changing.
SNAPSHOT_KEEP = None

_SCALARS = (bool, int, float, complex, str, bytes, type(None),
	np.bool_, np.integer, np.floating, np.complexfloating)

_MISS = object()

def _snapshot_value(val, containers = None):

	"""Return a snapshot of val, or _MISS if val should be shared not copied."""

	if containers is None:
		containers = SNAPSHOT_CONTAINERS

	if isinstance(val, np.ndarray):
		if val.size > SNAPSHOT_MAX_ELEMENTS:
			return _MISS
		return val.copy()

	if isinstance(val, _SCALARS):
		return val

	if isinstance(val, (list, tuple, dict, set, frozenset)):
		if containers == False:
			return _MISS
		try:
			return copy.deepcopy(val)
		except Exception:
			return _MISS

	return _MISS


def _class_state_keys(cls, use_skip = True):

	"""Names on the class that hold snapshot-worthy state."""

	keys = []

	for key, val in list(cls.__dict__.items()):
		if key.startswith('__'):
			continue
		if use_skip and (key in SNAPSHOT_SKIP):
			continue
		if callable(val):
			continue
		if isinstance(val, (property, classmethod, staticmethod)):
			continue
		keys.append(key)

	return keys


def snapshot_object_state(obj, keep = None, containers = None, use_skip = True):

	"""
	Capture the mutable numeric state of obj, instance level and class level.

	Returns a dict to hand back to restore_object_state. Heavy objects
	(interpolators, table queries, file handles) are deliberately not
	captured, because a proposal does not mutate them.
	"""

	cls = type(obj)

	if keep is None:
		keep = SNAPSHOT_KEEP

	inst = {}
	for key, val in list(obj.__dict__.items()):
		if use_skip and (key in SNAPSHOT_SKIP):
			continue
		if keep is not None and key not in keep:
			continue
		snap = _snapshot_value(val, containers)
		if snap is not _MISS:
			inst[key] = snap

	klass = {}
	for key in _class_state_keys(cls, use_skip = use_skip):
		if keep is not None and key not in keep:
			continue
		snap = _snapshot_value(cls.__dict__[key], containers)
		if snap is not _MISS:
			klass[key] = snap

	return {'instance': inst, 'class': klass, 'cls': cls}


def restore_object_state(obj, state):

	"""
	Write a snapshot back onto obj.

	Instance attributes are rebound. Class attributes are written in place
	when the shape still matches, so that a class-level array never gets
	shadowed by an instance-level one, which would silently split pide's
	state in two.
	"""

	for key, val in state['instance'].items():
		if isinstance(val, np.ndarray):
			cur = obj.__dict__.get(key, None)
			if isinstance(cur, np.ndarray) and cur.shape == val.shape:
				cur[:] = val
			else:
				obj.__dict__[key] = val.copy()
		elif isinstance(val, (list, tuple, dict, set, frozenset)):
			obj.__dict__[key] = copy.deepcopy(val)
		else:
			obj.__dict__[key] = val

	cls = state['cls']

	for key, val in state['class'].items():
		cur = cls.__dict__.get(key, None)
		if isinstance(val, np.ndarray):
			if isinstance(cur, np.ndarray) and cur.shape == val.shape:
				cur[:] = val
			else:
				setattr(cls, key, val.copy())
		elif isinstance(val, (list, tuple, dict, set, frozenset)):
			setattr(cls, key, copy.deepcopy(val))
		else:
			setattr(cls, key, val)

	#An attribute that appeared on the instance after the snapshot was taken
	#would otherwise survive a rejection. Remove any such leftovers.
	for key in [k for k in obj.__dict__ if k not in state['instance']]:
		if key in SNAPSHOT_SKIP:
			continue
		if _snapshot_value(obj.__dict__[key]) is not _MISS:
			del obj.__dict__[key]


def _values_match(a, b):

	if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
		a = np.asarray(a)
		b = np.asarray(b)
		if a.shape != b.shape:
			return False
		if a.dtype.kind in 'fc' and b.dtype.kind in 'fc':
			return bool(np.array_equal(a, b, equal_nan = True))
		return bool(np.array_equal(a, b))
	try:
		return bool(a == b)
	except Exception:
		return repr(a) == repr(b)


def check_state_restored(obj, state, raise_on_fail = True):

	"""
	Confirm that obj has been returned to exactly the state the snapshot
	captured. Call this straight after restore_object_state while testing.

	Returns the list of attribute names that did not match, which is empty
	when the rollback is complete.
	"""

	now = snapshot_object_state(obj)
	bad = []

	for scope in ('instance', 'class'):
		for key, val in state[scope].items():
			if key not in now[scope]:
				bad.append('%s.%s (missing after restore)' % (scope, key))
			elif not _values_match(val, now[scope][key]):
				bad.append('%s.%s' % (scope, key))
		for key in now[scope]:
			if key not in state[scope]:
				bad.append('%s.%s (appeared after restore)' % (scope, key))

	if bad and raise_on_fail:
		raise AssertionError('state not fully restored: ' + ', '.join(bad))

	return bad


def snapshot_report(obj, top = 15):

	"""
	Print what a snapshot copies and what it shares. Run this once before
	relying on the snapshot, and check that nothing a proposal mutates is
	in the shared list.
	"""

	cls = type(obj)
	copied = []
	shared = []

	def _classify(scope, key, val):
		snap = _snapshot_value(val)
		if snap is _MISS:
			shared.append((scope, key, type(val).__name__,
				val.nbytes if isinstance(val, np.ndarray) else -1))
		else:
			nbytes = val.nbytes if isinstance(val, np.ndarray) else 0
			copied.append((scope, key, type(val).__name__, nbytes))

	for key, val in list(obj.__dict__.items()):
		_classify('instance', key, val)
	for key in _class_state_keys(cls):
		_classify('class', key, cls.__dict__[key])

	total = sum(c[3] for c in copied)

	print('snapshot copies %d attributes, %.1f kB of array data'
		% (len(copied), total / 1024.0))
	print('snapshot shares %d attributes (not copied)' % len(shared))
	print()

	copied.sort(key = lambda c: -c[3])
	print('largest copied:')
	for scope, key, tname, nbytes in copied[:top]:
		print('  %-8s %-32s %-18s %8.1f kB' % (scope, key, tname, nbytes / 1024.0))

	print()
	print('shared (verify none of these change per proposal):')
	for scope, key, tname, nbytes in shared[:top]:
		extra = ('%.1f kB' % (nbytes / 1024.0)) if nbytes >= 0 else ''
		print('  %-8s %-32s %-18s %s' % (scope, key, tname, extra))
	if len(shared) > top:
		print('  ... and %d more' % (len(shared) - top))


def snapshot_profile(obj, repeat = 20, top = 20):

	"""
	Time the snapshot of each attribute separately and print the worst.

	Use this when snapshot_object_state is slower than expected. Anything
	near the top of the list that a proposal never changes belongs in
	SNAPSHOT_SKIP.
	"""

	import time

	cls = type(obj)
	rows = []

	def _time_one(scope, key, val):
		t0 = time.perf_counter()
		for _ in range(repeat):
			_snapshot_value(val)
		dt = (time.perf_counter() - t0) / repeat
		rows.append((dt, scope, key, type(val).__name__,
			len(val) if hasattr(val, '__len__') else 0))

	for key, val in list(obj.__dict__.items()):
		if key not in SNAPSHOT_SKIP:
			_time_one('instance', key, val)
	for key in _class_state_keys(cls):
		if key not in SNAPSHOT_SKIP:
			_time_one('class', key, cls.__dict__[key])

	rows.sort(key = lambda r: -r[0])
	total = sum(r[0] for r in rows)

	t0 = time.perf_counter()
	for _ in range(repeat):
		snapshot_object_state(obj)
	whole = (time.perf_counter() - t0) / repeat

	print('whole snapshot        %9.3f ms' % (whole * 1e3))
	print('sum of per-attribute  %9.3f ms  over %d attributes'
		% (total * 1e3, len(rows)))
	print()
	print('%-9s %-34s %-14s %8s %10s' % ('scope', 'attribute', 'type', 'len', 'ms'))
	for dt, scope, key, tname, n in rows[:top]:
		print('%-9s %-34s %-14s %8d %10.4f' % (scope, key, tname, n, dt * 1e3))

	print()
	print('to skip the worst offenders:')
	print('  from mcmc_state import SNAPSHOT_SKIP')
	print('  SNAPSHOT_SKIP.update([%s])'
		% ', '.join("'%s'" % r[2] for r in rows[:5]))


class ChangeWatcher(object):

	"""
	Record which attributes a proposal actually changes.

	Run this for a few thousand real iterations with SNAPSHOT_KEEP = None
	and SNAPSHOT_CONTAINERS = True, so that nothing is excluded from the
	comparison. It then prints the exact set of attributes that ever
	differed, which is the whitelist you want, measured rather than
	remembered.

	In _solv_MCMC_column:

		watcher = ChangeWatcher()                  # before the loop

		state_backup = snapshot_object_state(object)   # already there
		...
		watcher.compare(object, state_backup)      # after the forward model,
		                                           # before accept/reject

		watcher.report()                           # after the loop
	"""

	def __init__(self):
		self.counts = {}
		self.n = 0

	def compare(self, obj, state):

		self.n += 1
		now = snapshot_object_state(obj)

		for scope in ('instance', 'class'):
			for key in set(state[scope]) | set(now[scope]):
				a = state[scope].get(key, _MISS)
				b = now[scope].get(key, _MISS)
				if (a is _MISS) or (b is _MISS) or (not _values_match(a, b)):
					self.counts[key] = self.counts.get(key, 0) + 1

	def report(self, top = 200):

		rows = sorted(self.counts.items(), key = lambda r: -r[1])

		print('ChangeWatcher: %d proposals observed, %d attributes ever changed'
			% (self.n, len(rows)))
		print()
		print('%-40s %10s %9s' % ('attribute', 'changed', 'of total'))
		for key, c in rows[:top]:
			print('%-40s %10d %8.1f%%' % (key, c, 100.0 * c / max(self.n, 1)))

		print()
		print('SNAPSHOT_KEEP = set([')
		for key, c in rows:
			print("\t'%s'," % key)
		print('])')
		print()
		print('Anything changed by fewer than ~1%% of proposals is rare, not')
		print('optional. Keep every name printed above.')


def snapshot_dump(obj, n_ref = None, path = None):

	"""
	List every attribute on the object and its class, grouped, so you can
	pick a whitelist by deleting lines instead of recalling names.

	n_ref : int or None
		Length of the depth axis (len(object.T)). Arrays of this length are
		grouped first, since those are the per-depth fields a proposal is
		most likely to touch.
	path : str or None
		Write the pasteable SNAPSHOT_KEEP block to this file as well.
	"""

	cls = type(obj)

	if n_ref is None:
		T = getattr(obj, 'T', None)
		n_ref = len(T) if T is not None else -1

	groups = {
		'per-depth arrays (length %s)' % n_ref: [],
		'other arrays': [],
		'scalars and flags': [],
		'containers': [],
		'shared, not copyable': [],
	}

	def _place(scope, key, val):
		row = (scope, key, type(val).__name__,
			len(val) if hasattr(val, '__len__') else 0)
		if isinstance(val, np.ndarray):
			if val.ndim >= 1 and val.shape[0] == n_ref:
				groups['per-depth arrays (length %s)' % n_ref].append(row)
			else:
				groups['other arrays'].append(row)
		elif isinstance(val, _SCALARS):
			groups['scalars and flags'].append(row)
		elif isinstance(val, (list, tuple, dict, set, frozenset)):
			groups['containers'].append(row)
		else:
			groups['shared, not copyable'].append(row)

	for key, val in sorted(obj.__dict__.items()):
		_place('instance', key, val)
	for key in sorted(_class_state_keys(cls)):
		_place('class', key, cls.__dict__[key])

	lines = []
	for gname, rows in groups.items():
		print('=' * 70)
		print('%s   (%d)' % (gname, len(rows)))
		print('=' * 70)
		for scope, key, tname, n in rows:
			print('  %-9s %-36s %-14s len=%d' % (scope, key, tname, n))
		if gname != 'shared, not copyable':
			lines.append('\t# --- %s ---' % gname)
			for scope, key, tname, n in rows:
				lines.append("\t'%s'," % key)
		print()

	block = 'SNAPSHOT_KEEP = set([\n' + '\n'.join(lines) + '\n])\n'

	print('=' * 70)
	print('pasteable block, delete the lines you do not need')
	print('=' * 70)
	print(block)

	if path is not None:
		with open(path, 'w') as f:
			f.write(block)
		print('written to %s' % path)

	return block


class AdaptiveSnapshot(object):

	"""
	Learn which attributes a proposal changes, then snapshot only those.

	Phase 1, up to `learn_until` iterations: takes a full snapshot and
	records every attribute that differed after the proposal was applied.

	Phase 2, from then on: snapshots a narrowed set and skips containers,
	which is where nearly all of the cost is.

	The narrowed set is NOT just what was observed changing. It is the union
	of three things:

		1. every attribute observed to change during phase 1
		2. every array whose leading axis matches the depth axis
		3. every scalar, bool and string attribute

	Points 2 and 3 cost microseconds and make the set near-complete by
	construction, so a rare branch that never fired during learning cannot
	silently leak. Containers are the only thing genuinely dropped, and on
	a pide object those are static CSV tables or lists of references to
	arrays already covered by point 2.

	Usage in _solv_MCMC_column, three edits:

		tracker = AdaptiveSnapshot(learn_until = burning // 2)   # before loop

		state_backup = tracker.snapshot(object, _)               # replaces
		                                                         # snapshot_object_state

		tracker.observe(object, state_backup)                    # after the
		                                                         # forward model

	restore_object_state is unchanged, it works off whatever the state dict
	holds.
	"""

	def __init__(self, learn_until, n_depth = None, verbose = True):

		self.learn_until = int(learn_until)
		self.n_depth = n_depth
		self.verbose = verbose

		self.counts = {}
		self.n_observed = 0
		self.keep = None
		self.locked = False

	def _build_keep(self, obj):

		cls = type(obj)

		n_depth = self.n_depth
		if n_depth is None:
			T = getattr(obj, 'T', None)
			n_depth = len(T) if T is not None else -1

		keep = set(self.counts)
		n_observed_names = len(keep)

		def _consider(key, val):
			if isinstance(val, np.ndarray):
				if val.ndim >= 1 and val.shape[0] == n_depth:
					keep.add(key)
			elif isinstance(val, _SCALARS):
				keep.add(key)

		for key, val in list(obj.__dict__.items()):
			if key not in SNAPSHOT_SKIP:
				_consider(key, val)
		for key in _class_state_keys(cls):
			_consider(key, cls.__dict__[key])

		if self.verbose:
			print('AdaptiveSnapshot: locked after %d observed proposals' % self.n_observed)
			print('  %d attributes changed during learning' % n_observed_names)
			print('  %d attributes kept (changed + per-depth arrays + scalars)' % len(keep))
			rare = sorted([(c, k) for k, c in self.counts.items() if c <= max(1, self.n_observed // 100)])
			if len(rare) > 0:
				print('  rarely changed, kept: %s'
					% ', '.join('%s(%d)' % (k, c) for c, k in rare[:10]))

		return keep

	def snapshot(self, obj, iteration):

		if (self.locked == False) and (iteration >= self.learn_until):
			self.keep = self._build_keep(obj)
			self.locked = True

		if self.locked == True:
			return snapshot_object_state(obj, keep = self.keep, containers = False)

		return snapshot_object_state(obj, keep = None, containers = True)

	def observe(self, obj, state):

		if self.locked == True:
			return

		self.n_observed += 1
		now = snapshot_object_state(obj, keep = None, containers = True)

		for scope in ('instance', 'class'):
			for key in set(state[scope]) | set(now[scope]):
				a = state[scope].get(key, _MISS)
				b = now[scope].get(key, _MISS)
				if (a is _MISS) or (b is _MISS) or (not _values_match(a, b)):
					self.counts[key] = self.counts.get(key, 0) + 1


def snapshot_cost_report(obj, changed = None, top = 40, cutoff = 0.98, repeat = 20):

	"""
	Rank every attribute by how long it takes to snapshot, so you can skip
	the expensive ones.

	changed : iterable of str or None
		Names known to change during a proposal, e.g. the list ChangeWatcher
		printed. These are marked KEEP and are never suggested for skipping,
		however expensive they are. Pass it. Without it the suggested skip
		list is ranked purely by cost and will happily propose dropping real
		state.
	top : int
		How many rows to print.
	cutoff : float
		Build the suggested skip list from the costliest attributes until
		this fraction of the skippable cost is covered.

	Prints a table and a pasteable SNAPSHOT_SKIP block.
	"""

	import time

	cls = type(obj)
	changed = set(changed) if changed is not None else set()

	rows = []

	def _time_one(scope, key, val):
		t0 = time.perf_counter()
		for _ in range(repeat):
			_snapshot_value(val, containers = True)
		dt = (time.perf_counter() - t0) / repeat
		rows.append([dt, scope, key, type(val).__name__,
			len(val) if hasattr(val, '__len__') else 0])

	for key, val in list(obj.__dict__.items()):
		if key not in SNAPSHOT_SKIP:
			_time_one('instance', key, val)
	for key in _class_state_keys(cls):
		if key not in SNAPSHOT_SKIP:
			_time_one('class', key, cls.__dict__[key])

	rows.sort(key = lambda r: -r[0])
	total = sum(r[0] for r in rows)
	keep_cost = sum(r[0] for r in rows if r[2] in changed)

	print('%d attributes, %.3f ms total' % (len(rows), total * 1e3))
	if len(changed) > 0:
		print('%.3f ms of that is in attributes that change (must be kept)'
			% (keep_cost * 1e3))
		print('%.3f ms is skippable' % ((total - keep_cost) * 1e3))
	print()
	print('%4s %-9s %-34s %-12s %8s %9s %8s %6s'
		% ('#', 'scope', 'attribute', 'type', 'len', 'ms', 'cum%', 'note'))

	run = 0.0
	for i, (dt, scope, key, tname, n) in enumerate(rows[:top]):
		run += dt
		note = 'KEEP' if key in changed else ''
		print('%4d %-9s %-34s %-12s %8d %9.4f %8.1f %6s'
			% (i + 1, scope, key, tname, n, dt * 1e3, 100.0 * run / total, note))

	skippable = [r for r in rows if r[2] not in changed]
	budget = (total - keep_cost) * cutoff
	chosen = []
	run = 0.0
	for dt, scope, key, tname, n in skippable:
		if run >= budget:
			break
		chosen.append(key)
		run += dt

	print()
	print('skipping the %d below removes %.3f ms of %.3f ms skippable'
		% (len(chosen), run * 1e3, (total - keep_cost) * 1e3))
	print()
	print('SNAPSHOT_SKIP = set([')
	line = []
	for key in sorted(chosen):
		line.append("'%s'," % key)
		if len(line) == 4:
			print('\t' + ' '.join(line))
			line = []
	if len(line) > 0:
		print('\t' + ' '.join(line))
	print('])')

	if len(changed) == 0:
		print()
		print('WARNING: no `changed` set was given, so this list is ranked by')
		print('cost alone and may contain attributes a proposal modifies.')
		print('Run ChangeWatcher first and pass its names as `changed`.')

	return chosen


class RollbackAudit(object):

	"""
	End-to-end test that a rejected proposal leaves the object exactly as it
	was, including the attributes SNAPSHOT_SKIP excludes.

	check_state_restored compares against the snapshot, so it cannot see a
	name the snapshot never captured. This does the opposite: it takes a
	FULL picture before the proposal, ignoring SNAPSHOT_SKIP entirely, and
	another after the restore, and reports anything that differs. That is
	the only check that can catch a wrongly skipped attribute.

	Three edits, all temporary:

		audit = RollbackAudit()                      # before the loop

		audit.before(object)                         # next to
		                                             # state_backup = ...

		restore_object_state(object, state_backup)   # already there
		audit.after(object)                          # in BOTH rejection
		                                             # branches

		audit.report()                               # after the loop

	Expect it to be slow, roughly 30 ms per rejected proposal, since it
	takes two unfiltered snapshots. Run it for a couple of thousand
	iterations, read the report, then take it out again.
	"""

	def __init__(self):
		self._pre = None
		self.counts = {}
		self.n_checked = 0

	def before(self, obj):
		self._pre = snapshot_object_state(obj, keep = None,
			containers = True, use_skip = False)

	def after(self, obj):

		if self._pre is None:
			return

		self.n_checked += 1
		post = snapshot_object_state(obj, keep = None,
			containers = True, use_skip = False)

		for scope in ('instance', 'class'):
			for key in set(self._pre[scope]) | set(post[scope]):
				a = self._pre[scope].get(key, _MISS)
				b = post[scope].get(key, _MISS)
				if (a is _MISS) or (b is _MISS) or (not _values_match(a, b)):
					self.counts[key] = self.counts.get(key, 0) + 1

		self._pre = None

	def report(self):

		print('RollbackAudit: %d rejections checked' % self.n_checked)

		if self.n_checked == 0:
			print('  nothing was checked. after() never ran, so no proposal')
			print('  was rejected, or the call is in the wrong branch.')
			return

		if len(self.counts) == 0:
			print('  PASS, every rejection restored the object exactly,')
			print('  including everything in SNAPSHOT_SKIP.')
			return

		print('  FAIL, %d attributes did not come back:' % len(self.counts))
		print()
		print('  %-38s %10s %9s %s' % ('attribute', 'leaked', 'of checked', 'in SNAPSHOT_SKIP'))
		for key, c in sorted(self.counts.items(), key = lambda r: -r[1]):
			print('  %-38s %10d %8.1f%%  %s'
				% (key, c, 100.0 * c / self.n_checked,
				'YES -> remove it from the set' if key in SNAPSHOT_SKIP else 'no'))
		print()
		print('  A name marked YES is being skipped but does change and is not')
		print('  an alias of an array that is captured. Take it out of')
		print('  SNAPSHOT_SKIP. A name marked no means restore_object_state')
		print('  itself failed on it, which is a different problem.')