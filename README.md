# icon-pamtra-insitu
This is an open copy of pamtra-insitu codes from the WarmWorld-Easier project "insitu"

# Catalog for Pamtra insitu scripts with YAC

# Setup on levante
To setup the pamtra-insitu workflow on levante follow the following steps:

Build ICON
```bash
curl -s -L  https://gitlab.dkrz.de/icon/icon-model/-/archive/release-2024.07-public/icon-model-release-2024.07-public.tar.gz | tar xzf -
mv icon-model-release-2024.07-public/ icon
cd icon
mkdir build
cd build
../config/dkrz/levante.gcc --enable-opnemp --enable-python-bindings --with-pic
make -j32
pip install externals/yac/python
```

# Prepare runscripts
To attach a pamtra component to ICON you need to modify the runscript as follows:

To enable the "output_coupling" add `coupled_to_model = .TRUE.` to the `coupling_mode_nml`:
```
&coupling_mode_nml
  coupled_to_model = .TRUE.
/
```

The we need to add an additional python process to MPI. We use the `--multi-prog` feature of `srun`.
Add this just before the `srun` (`${START}`) call
```bash
cat > mpmd.conf << EOF
/path/to/venv/bin/python path/to/pamtra/aircraft_nopamtra.py --track=TRACKFILE.nc
* ${MODEL}
EOF
```
(adjust the pathes) and replace the line
```
${START} ${MODEL}
```
with
```
${START} --multi-prog mpmd.conf
```
(This does only work with slurm.)

## TRACKFILE.nc
This is just a netCDF file that contains the time,lat,lon information of the satellite track. This is the ncdump of my file:
```
dimensions:
	time = 17523 ;
variables:
	int time(time) ;
		time:units = "seconds since 2017-01-01T00:00:00+00:00" ;
		time:calendar = "proleptic_gregorian" ;
	double seconds(time) ;
		seconds:_FillValue = NaN ;
		seconds:long_name = "Seconds since start of day" ;
		seconds:units = "s" ;
		seconds:comment = "Time as recorded by the GPS" ;
	double lat(time) ;
		lat:_FillValue = NaN ;
		lat:long_name = "latitude" ;
		lat:units = "degrees_north" ;
		lat:comment = "GPS latitude" ;
	double lon(time) ;
		lon:_FillValue = NaN ;
		lon:long_name = "longitude" ;
		lon:units = "degrees_east" ;
		lon:comment = "GPS longitude" ;
	double alt(time) ;
		alt:_FillValue = NaN ;
		alt:long_name = "altitude" ;
		alt:units = "m" ;
		alt:comment = "Aircraft altitude" ;
```
but the variables "seconds" and "alt" are not relevant
