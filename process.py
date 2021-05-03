from sr.sim_cfg import SimpleCfg
from sr.volvox import analyze, load_data
from sr.volvox.io import get_xyz_data
from sr.volvox.process import get_orientation_vec
from pathlib import Path
import json
import pandas as pd

import matplotlib.pylab as plt
from sr.numpyc import *
sims =  set([sim.parent for sim in Path(".").glob("**/*.san")])
results = {}
results_old = pd.read_hdf('results2.h5', 'v1')

for idx, wdir in enumerate(sims):
    print(f'procesing {wdir.name} ({idx}/{len(sims)}) {idx/len(sims)*100:.0f}')
    cfg = SimpleCfg(wdir=wdir)
    cfg.load_from_file()

    if str(wdir) in results_old.index:
        if results_old.loc[str(wdir)].globalt >= cfg._cfg.globalt:
            print('==== skipping ====')
            continue
        else:
            print('UPDATE')
    try:
        with open(wdir.joinpath("parameters.json")) as fp:
            params = json.load(fp)
        data = load_data(wdir)
    except Exception as e:
        print(e)
        continue
    
    try:
        res, N, V = analyze(data)
    except Exception as e:
        print(e)
        continue
    params.update(res)
    params['globalt'] = data.last_valid_index()

    
    if cfg._cfg.globalt - data.last_valid_index() > 250.0:
        print('ERROR with {wdir}')
    params['tau_power'] = cfg._cfg.md_objects[0].cilia.power_stroke_duration 
    params['tau_recover'] = cfg._cfg.md_objects[0].cilia.recovery_stroke_duration  
    params['tau'] = params['tau_power'] + params['tau_recover']
    params['num_cilia'] = sum(cfg._cfg.md_objects[0].num_equators)
    params['thetaZ'] = cfg._cfg.md_objects[0].cilia.thetaZ
    params['run_type'] = wdir.parts[0]
    results[str(wdir)] = params

if results:
    results = pd.DataFrame(results).T
    results["alpha"] = results.a/np.pi*180
    results.to_hdf('results2.h5', 'v1b')
