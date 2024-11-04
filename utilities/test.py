import _ctypes
# test_ctypes.py
try:
    import _ctypes
    print("_ctypes is available")
except ModuleNotFoundError:
    print("_ctypes is not available")

