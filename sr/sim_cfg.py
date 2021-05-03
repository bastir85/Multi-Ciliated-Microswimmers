from pathlib import PurePath

import pylibconfig2 as cfg
from tables import HDF5ExtError, NoSuchNodeError
import numpy

class MDObject(object):

    def __init__(self, name):
        self.name = name
        self._data = None
        self._cfg = None

    def _get_data(self, attr):
        attr = "_" + attr
        class_attr = self.__dict__.get(attr)
        if class_attr is None:
            #print("No data loaded. Use load_h5!!")
            return None
        return class_attr

    def __ne__(self, other):
        return self.name != other.name

    def __eq__(self, other):
        return self.name == other.name
    
    def __hash__(self):
        return hash(self.name)
        
    def __getattribute__(self, attr):
        cfg = object.__getattribute__(self, "_cfg")
        if cfg:
            cfg_attr = cfg.get(attr)
            if cfg_attr is not None:
                return cfg_attr

        #default behavior
        try:
            class_attr = object.__getattribute__(self, attr)
            return class_attr
        except AttributeError:
            pass
        return None
        #data_attr = self._get_data(attr)
        #if data_attr is None:
        #    pass
        #    #raise ValueError
       
        #return data_attr
        


    def to_dict(self, force=False):
        '''
            Be carful this overwrites the original if force is true
        '''
        if self._cfg and not force:
            return self._cfg

        parameters = dict([(k,v) for k,v in self.__dict__.items() if not k.startswith("_") ])
        if self.position:
            parameters['position'] = cfg.ConfArray(self.position)
        if self.y_thetas:
            parameters['y_thetas'] = cfg.ConfArray(self.y_thetas)
    
        if self._cfg:
            self._cfg.__dict__.update(parameters)
        else:
            self._cfg = cfg.ConfGroup(parameters)
        
        return self._cfg

    def from_cfg(self, c_obj):
        self._cfg = c_obj


    @property
    def data(self):
        return self._data


class BasicCfg(object):
    def __init__(self, wdir): # Used to take everything.
        self.md_objects= set() 
        self.parameters={}
        self.wdir = wdir

    def add(self, obj):
        self.md_objects.add(obj)

    def load_from_file(self, filename="sim.cfg"): 
        fp = open(PurePath(self.wdir,filename))
        cfg_str = "\n".join(fp.readlines())
        fp.close()
        self._load_cfg(cfg_str)

    def _load_cfg(self, cfg_str):
        config = cfg.Config(cfg_str)
        self.md_objects = set()
        for c_obj in config.md_objects:
            obj = MDObject(c_obj.name)
            obj.from_cfg(c_obj)
            self.md_objects.add(obj)
        self.md_objects = list(self.md_objects)
        self._cfg = config 
        
    def _update_fields(self):    
        for k, v in self._cfg.__dict__.items():
            if k == "md_objects": continue
            self.parameters[k] = v 
            self.__dict__[k] =v 
    
    @property
    def data(self):
        if hasattr(self,"_data"):
            return self._data
        if len(self.md_objects) == 1:
            return self.md_objects[0]._data
    @property
    def md_obj(self):
        return self.md_objects[0]

class SimpleCfg(BasicCfg):

    def save_cfg(self, filename):
        strdata = self._get_config()
        fp = open(filename, 'w')
        fp.write(str(strdata))
        fp.close()

    def _get_config(self, force=False):
        config = cfg.Config()
        for k,v in self.parameters.items():
            if type(v) == numpy.float64:
                v = float(v)
            config.set(k,v)
        c_md_objects = cfg.ConfList([obj.to_dict(force=force) for obj in self.md_objects])
        config.setup("md_objects", cfg.ConfList(c_md_objects))
        return config

