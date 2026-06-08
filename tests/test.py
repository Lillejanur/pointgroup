import sys

print(sys.executable)

try:
    import ase
    print("ASE:", ase.__version__)
except Exception as e:
    print("ASE error:", e)