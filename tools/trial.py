#!/usr/bin/env python3

import os

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../seel_src')
print(core_path_ext)

from cond_models.melt_odd import * 
from cond_models.fluids_odd import * 
#importing odd rock functions
from cond_models.rocks.granite_odd import * 
from cond_models.rocks.granulite_odd import *
from cond_models.rocks.sandstone_odd import *
from cond_models.rocks.gneiss_odd import *
from cond_models.rocks.amphibolite_odd import *
from cond_models.rocks.basalt_odd import *
from cond_models.rocks.mud_odd import *
from cond_models.rocks.gabbro_odd import *
from cond_models.rocks.other_rocks_odd import *
#importing odd mineral functions
from cond_models.minerals.quartz_odd import *
from cond_models.minerals.plag_odd import *
from cond_models.minerals.amp_odd import *
from cond_models.minerals.kfelds_odd import *
from cond_models.minerals.opx_odd import *
from cond_models.minerals.cpx_odd import *
from cond_models.minerals.mica_odd import *
from cond_models.minerals.garnet_odd import *
from cond_models.minerals.ol_odd import *
from cond_models.minerals.mixtures_odd import *
from cond_models.minerals.other_odd import *