# Make turbulent velocity field
Captured from https://github.com/mikegrudic/MakeCloud .

## How to run:

Determin `res` in `make_turb.py` script (e.g. 64) which makes three velocity cubes of res^3 file. Then run:

```sh
python ./make_turb.py
gfortran ./convert_python_cube.f90 -o ./convert_python_cube.exe
./convert_python_cube.exe
```

This will make three velocity cubes that can be used as `phantom` input. Just copy them in `phantom/data/velfield` directory.
