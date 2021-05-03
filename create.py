# coding: utf-8

from sr.sim_cfg import  *
from math import pi, ceil
import random
import string
import os
import numpy as np
import pandas as pd
import pickle
import pylibconfig2 as cfg2
import json

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))


class Cilia(MDObject):
    def __init__(self, max_energy, start_delay=0.0, thetaZ=0.0):
        springconstant = 20000.0
        self.bending = float(np.sqrt(2*max_energy/springconstant))
        self.springconstant = springconstant 
        self.beatstrength = 0.17;
        self.angle = 0.0;
        self.name = "cilia";
        self.rodlength = 26;
        self.thetaZ = float(thetaZ);
        self.thetaX = 0.0;
        self.thetaY = 0.0;
        self.m_speed = 5.3;
        self.max_bend = 0.7;
        self.position = [ 0.0, 0.0, 0.0];
        self.min_bend = -0.7;
        self.write_data = False;
        
        
        self.power_stroke_duration =  75.0;
        self.recovery_stroke_duration = 150.0;
        self.start_delay = start_delay
        super().__init__("ConstantBeatCilia")
        self.type = "ConstantBeatCilia"
        self._taub = self.power_stroke_duration + self.recovery_stroke_duration 


def get_delays(phi_rings, num_cilia_per_ring, k_theta, theta_offsets, phase_lag):
    assert k_theta == int(k_theta)
    R = 8
    ## k_theta = np.sin(a)
    k_phi = phase_lag*4/np.pi
    res = []
    for idx, phi in enumerate(phi_rings):
        for theta in np.arange(num_cilia_per_ring[idx])/num_cilia_per_ring[idx]*2*np.pi+ theta_offsets[idx]:
            res.append((k_phi*phi+k_theta*theta))
    return np.array(res) % (2*np.pi)


def start_sim(proj_dir, k_theta, energy=5,  equator_num=3, thetaZ=0.0, phase_lack=None):
    ring2_theta = 45.0
    ring2_num = int(np.round(equator_num*np.cos(ring2_theta/180*np.pi)))
    if ring2_num == 0:
        raise ValueError(f"Density on ring {ring2_theta} is zero, abort!")

    name = f"volvox_E{energy:.0f}_ktheta{k_theta:.3f}"
    name +='_'+''.join(random.choice(string.ascii_letters + string.digits) 
                                            for _ in range(2))
    wdir = os.path.join(ROOT_DIR, proj_dir, f"n_eq{equator_num}", f"thetaZ{thetaZ/np.pi*180}" ,f"k_theta{k_theta}", name)
    os.makedirs(wdir)
    os.chdir(wdir)

    sim = SimpleCfg(name) 
    sim.parameters['name'] = "volvox"
    sim.parameters['globalt'] = 0.
    sim.parameters['tstep'] = 0.05 
    sim.parameters['runtime'] = 23.8*60**2
    sim.parameters['sizex'] = 100
    sim.parameters['sizey'] = 100
    sim.parameters['sizez'] = 100
    sim.parameters['omega'] = 0.0
    
    volvox = MDObject('volvox')
    cilia = Cilia(max_energy=energy, thetaZ=thetaZ)
    volvox.cilia = cilia.to_dict()
    volvox.position = [0., 0., 0.]
    volvox.type = "Clami"
    volvox.springconstant = 20000.0
    volvox.radius = 8.0
    volvox.y_thetas = cfg2.ConfArray([-ring2_theta, 0., ring2_theta])
    num_cilia_per_ring = [ring2_num, equator_num, ring2_num]
    volvox.num_equators = cfg2.ConfArray(num_cilia_per_ring)
    delta = 2*np.pi/equator_num/2
    volvox.theta_offset = cfg2.ConfArray([delta, 0.0, delta])
    if True:
        delays = np.round(
            get_delays(phi_rings=np.array(volvox.y_thetas)/180*np.pi, num_cilia_per_ring=num_cilia_per_ring,
                       k_theta=k_theta,
                       phase_lag=phase_lack,
                       theta_offsets=list(volvox.theta_offset))/2/np.pi*cilia._taub, 2)
    else:
        delays = np.zeros(np.sum(num_cilia_per_ring))
    volvox.start_delays = cfg2.ConfArray([float(d) for d in delays])
    volvox.cilia.min_bend = -0.7 
    volvox.cilia.max_bend = 0.7 
    sim.add(volvox)
    
    sim.parameters['desiredt'] = float(ceil(1000*cilia._taub))
 
    sim.save_cfg("sim.cfg")
    with open('parameters.json', 'w') as fp:
        k_phi = phase_lack*4/np.pi
        json.dump({'a': np.arctan2(k_theta, k_phi), 'phase_lack': phase_lack, 'k_theta': k_theta,
                   'thetaZ': thetaZ, 'num_cil_eq': equator_num,
                   'num_cilia': int(np.sum(num_cilia_per_ring)),
                   'max_energy': energy}, fp)


wdir="final_long_smooth"
q = np.linspace(-np.pi, np.pi, 7,endpoint=False)
phase_lacks = q + 2*np.pi/14 

idx = 0
num_cil_eq = 10 

if __name__ == "__main__":
    start_sim(wdir, k_theta=1, thetaZ=0.0, energy=1.0, equator_num=num_cil_eq, phase_lack=phase_lacks[0])

def test():
    for phase_lack in phase_lacks:
        for k_theta in range(0, num_cil_eq//2+1):
            start_sim(wdir, k_theta=k_theta, thetaZ=0.0, energy=1.0, equator_num=num_cil_eq, phase_lack=phase_lack)
        
        for thetaZ in np.linspace(0, np.pi/4, 5+2)[1:]: 
            start_sim(wdir, k_theta=0.0, thetaZ=thetaZ, energy=1.0, equator_num=num_cil_eq, phase_lack=phase_lack)

