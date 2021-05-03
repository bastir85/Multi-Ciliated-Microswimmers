import numpy as np
import pandas as pd

def _iter_loadxyz(infile, skiprows=0, dtype=float):
    def iter_func():
        for _ in range(skiprows):
            next(infile)
        frame = None
        globalT = None
        for line in infile:
            line = line.strip().split()
            if len(line) == 1:
                try:
                    num_parts = int(line[0])
                except ValueError:
                    print("H")
                    print(line)


                if frame:
                    if len(frame) == num_parts:
                        yield globalT, frame
                    else:
                        print ("Frame {0} defect.".format(globalT))
                        yield globalT, [[None]*3]*num_parts
                frame = []
                try:
                    line = next(infile)
                    globalT = float(line.strip().split()[-1])
                except (IndexError, ValueError):
                    if globalT is None:
                        globalT = 0
                    else:
                        globalT += 1 ##None
                    #print("T")
                    #print(line)
                continue
            vec = []
            for item in line[1:]:
                vec.append(dtype(item))
            frame.append(vec)
        yield globalT, frame
            #iter_loadxyz.rowlength = len(line)
    return iter_func()


def load_xyz(filename):
    import numpy as np
    import pandas as pd
    data_xyz = _iter_loadxyz(filename)
    tsteps, frames = zip(*data_xyz)
    frames = np.array(frames)
    frames = frames.swapaxes(1,2).reshape(-1,frames.shape[1])
    idx = pd.MultiIndex.from_product([tsteps,[0,1,2]], names=['time', 'coordinate'])
    data = pd.DataFrame(frames, index=idx)
    return data


def load_data(wdir):
    columns = ["time", "x", "y", "z", "vx", "vy", "vz",
               "ekin",
               "nx", "ny", "nz",
               "bx", "by", "bz",
               "px", "py", "pz",
               "power"]
    data =  pd.read_csv(wdir.joinpath("volvox.san"), sep=" ", 
                        header=None, names=columns)
    data.index= data.time#pd.to_timedelta(data.time, unit='s')
    ncols = data.dropna(axis=1).shape[1]

    # v1 version
    if ncols == 11:
        data.dropna(inplace=True, axis=1)
        data.columns =  ['time', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'nx', 'ny', 'nz',
                         'power']
    elif all(data.power.isna()):
        data.dropna(inplace=True, axis=1)
        data.columns =  ['time', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'nx', 'ny', 'nz', 'bx',
               'by', 'bz', 'px', 'py', 'pz', 'power']

    for col_name in set(columns) - set(data.columns):
        data[col_name] = np.nan

    return data


def get_xyz_data(wdir):
    with open(wdir.joinpath('volvox.xyz')) as fp:
        data = load_xyz(fp)
    sphere = data.loc[:,:642].unstack().stack(level=0)
    cilias = data.loc[:,643:].unstack().stack(level=0)
    cilias['cilidx'] =  cilias.groupby('time').apply(lambda l: l.groupby(np.arange(len(l))//26).ngroup()).droplevel(0)
    cilias = cilias.reset_index().set_index(['time','cilidx','level_1'], drop=True)
    cilias = cilias.loc[pd.IndexSlice[:,::3,:], : ].reset_index()
    cilias['cilidx'] = cilias.cilidx //3
    cilias = cilias.set_index(['time','cilidx','level_1'], drop=True)
    return sphere, cilias


def get_xyz_data_reduced(wdir, reduced_rodlength=39):
    with open(wdir.joinpath('volvox.xyz')) as fp:
        data = load_xyz(fp)
    cilias = data.unstack().stack(level=0)
    cilias['cilidx'] =  cilias.groupby('time').apply(lambda l: l.groupby(np.arange(len(l))//reduced_rodlength).ngroup()).droplevel(0)
    cilias = cilias.reset_index().set_index(['time','cilidx','level_1'], drop=True)
    cilias.index.names =  ['time', 'cilidx', 's']
    return cilias

