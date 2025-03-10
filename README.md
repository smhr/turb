# Make turbulent velocity field
Captured from https://github.com/mikegrudic/MakeCloud .

## How to run:

Determin `res` in `make_turb.py` script (e.g. 64). Then run:

```sh
python ./make_turb.py
```

This makes three velocity cubes of res^3 file: `vx.bin, vy.bin, vz.bin`.

Then convert them to an appropriate format for `phantom`:

```sh
gfortran ./convert_python_cube.f90 -o ./convert_python_cube.exe
./convert_python_cube.exe
```

This will make three velocity cubes `cube_v1.dat, cube_v2.dat, cube_v3.dat` that can be used as `phantom` input. Just copy them in `phantom/data/velfield` directory.

## How to plot turbulent velocity field

```sh
cd test
python3 ./plot_v_field.py
```