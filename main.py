from topopt import main
from dolfin import *
set_log_level(50)


# The real main driver
if __name__ == "__main__":
    main(nelx=60, nely=20, volfrac=0.5, penal=3.0, rmin=2.0)
