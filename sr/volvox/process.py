import numpy as np

def get_rotation_vector_ts_p(e,b,p, time_step,k=1):
    import pandas as pd
    We = list()
    Wb = list()
    Wp = list()
    r = -(len(e) % k) or None
    for l in range(k):
        Ie = np.cross(e.iloc[l:r:k], np.gradient(e.iloc[l:r:k], time_step*k, axis=0))
        Ib = np.cross(b.iloc[l:r:k], np.gradient(b.iloc[l:r:k], time_step*k, axis=0))
        Ip = np.cross(p.iloc[l:r:k], np.gradient(p.iloc[l:r:k], time_step*k, axis=0))
        We.append(np.sum(Ip*e.iloc[l:r:k],axis=1))
        Wb.append(np.sum(Ip*b.iloc[l:r:k],axis=1))
        Wp.append(np.sum(Ie*p.iloc[l:r:k],axis=1))
    We = pd.concat(We).sort_index()
    Wb = pd.concat(Wb).sort_index()
    Wp = pd.concat(Wp).sort_index()
    return pd.DataFrame({"Omega_n": We, "Omega_b":Wb, "Omega_p":Wp})


def normalize(vec):
    return vec.divide(np.sqrt((vec**2).sum(axis=1)), axis=0)

def get_orientation_vec(data):
    N = normalize(data[["nx", "ny", "nz"]])
    B = normalize(data[["bx", "by", "bz"]])
    P = normalize(data[["px", "py", "pz"]])
    return N, B, P


def analyze(data, eq_time=500.0):
    tstep = data.time.diff().mean() 
    assert tstep == 2.5

    results = {}
    data = data.loc[eq_time:]
    N, B, P = get_orientation_vec(data)
    phi = np.arccos((N*N.iloc[0]).sum(axis=1))
    results['phi_mean'] = phi.mean()
    results['phi_std'] = phi.std()
    results['phi_max'] = phi.max()
    
    V = data[["vx", "vy", "vz"]]
    results['vx_m'] = V.vx.mean()
    results['vy_m'] = V.vy.mean()
    results['vz_m'] = V.vz.mean()
    
    results['v_abs_x_m'] = V.vx.abs().mean()
    results['v_abs_y_m'] = V.vx.abs().mean()
    results['v_abs_z_m'] = V.vx.abs().mean()
    
    results['v_abs_x_std'] = V.vx.abs().std()
    results['v_abs_y_std'] = V.vx.abs().std()
    results['v_abs_z_std'] = V.vx.abs().std()
    
    results['v_normal_m'] = (N*V.values).sum(axis=1).mean()
    results['v_normal_std'] = (N*V.values).sum(axis=1).std()

    Omega = get_rotation_vector_ts_p(N,B,P, tstep, 240)
    results['Omega_n_mean'] = Omega.Omega_n.mean()
    results['Omega_n_std'] = Omega.Omega_n.std()
    results['Omega_b_mean'] = Omega.Omega_b.mean()
    results['Omega_b_std'] = Omega.Omega_b.std()
    results['Omega_p_mean'] = Omega.Omega_p.mean()
    results['Omega_p_std'] = Omega.Omega_p.std()

    E_kin = np.sqrt((V*V).sum(axis=1))
    results['e_kin_m'] = E_kin.mean()
    results['e_kin_std'] = E_kin.std()
    
    results['power_total'] = data.power.sum()
    results['power_mean'] = data.power.mean()
    results['power_std'] = data.power.std()
    
    return results, N, V
