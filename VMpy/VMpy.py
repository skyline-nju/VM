# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

"""
 Simulating Vicsek model in 2D.

    Caution: when calling get_coarse_grained_snap(num, svx, svy, l) in python,
    the data type of num should be np.int32 otherwise a error will occur in Linux
    platform.

"""


from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_VMpy')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_VMpy')
    _VMpy = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_VMpy', [dirname(__file__)])
        except ImportError:
            import _VMpy
            return _VMpy
        try:
            _mod = imp.load_module('_VMpy', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _VMpy = swig_import_helper()
    del swig_import_helper
else:
    import _VMpy
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0


def set_random_seed(seed):
    """set_random_seed(int seed)"""
    return _VMpy.set_random_seed(seed)

def set_v0(_v0):
    """set_v0(double _v0)"""
    return _VMpy.set_v0(_v0)

def set_eta(_eta):
    """set_eta(double _eta)"""
    return _VMpy.set_eta(_eta)

def set_eps(_eps):
    """set_eps(double _eps)"""
    return _VMpy.set_eps(_eps)

def setLx(_Lx):
    """setLx(double _Lx)"""
    return _VMpy.setLx(_Lx)

def setLy(_Ly):
    """setLy(double _Ly)"""
    return _VMpy.setLy(_Ly)

def ini(*args):
    """
    ini(double * x, double * y, double * vx, double * vy, int seed, double v0, double _eta, double _eps, double Lx, double Ly)
    ini(double * x, double * y, double * vx, double * vy, int nBird, int seed, double v0, double _eta, double _eps, double Lx, double Ly)
    """
    return _VMpy.ini(*args)

def run(nstep):
    """run(int nstep)"""
    return _VMpy.run(nstep)

def get_snap(*args):
    """
    get_snap(double * x, double * y, double * vx, double * vy)
    get_snap(double * x, double * y, double * vx, double * vy, int nBird)
    """
    return _VMpy.get_snap(*args)

def get_coarse_grained_snap(*args):
    """
    get_coarse_grained_snap(int * num, double * svx, double * svy, double l)
    get_coarse_grained_snap(int * num, double * svx, double * svy, int ncells, double l)
    """
    return _VMpy.get_coarse_grained_snap(*args)
# This file is compatible with both classic and new-style classes.


