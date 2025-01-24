import numpy as np
import xarray as xr
import pandas as pd
import sys
from argparse import ArgumentParser, FileType
from pathlib import Path

this_dir = Path(__file__).parent

parser = ArgumentParser()
parser.add_argument("--descriptorfile", type=FileType('r'), default=this_dir.parent/"descriptorfiles"/"icon_1mom.py", help="")
parser.add_argument("--track", type=str, help="File containing time-lon-lat data") # TODO Height?
parser.add_argument("--frequency", type=float, nargs="+", default=94)
parser.add_argument("--dt", type=str, default=None, help="Timestep") #TODO output time? for now I stay strict to model timestep and assume mean aircraft location ...
parser.add_argument("--source_component", type=str, default="atmo_output", help="component name of the source component")
parser.add_argument("--source_grid", type=str, default="icon_atmos_grid", help="grid name of the source grid")
args = parser.parse_args()


from yac import YAC, UnstructuredGrid, Location, Field, TimeUnit, \
    InterpolationStack, NNNReductionType, Reduction, Action
import pyPamtra


yac = YAC()

comp = yac.def_comp("pamtra")

track = xr.open_dataset(args.track)
lats = np.deg2rad(track.lat.values)
lons = np.deg2rad(track.lon.values) # TODO knowing the model timesteps I can in principle calculate the mean locations in advance and form a reduced grid, maybe pass as an argument?
                                    # Also, knowing the model domain I can exclude the track points outside the domain in advance?? probably not worth the complication
# alts = track.alt.values # need to think about this
times = track.time.values
ncoords = len(times)

grid = UnstructuredGrid("pamtra-grid", [3]*ncoords,
                        [*lons, *(lons+0.1), *lons],
                        [*lats, *lats, *(lats+0.1)],
                        np.array(list(zip(range(ncoords), range(ncoords, 2*ncoords), range(2*ncoords, 3*ncoords)))).ravel())
points = grid.def_points(Location.CELL,
                         lons, lats)

yac.sync_def()

# print(yac.get_field_names(*source), file=sys.stderr)

interp = InterpolationStack()
interp.add_nnn(NNNReductionType.AVG, 1, 1.0)

source = (args.source_component, args.source_grid)

dt = yac.get_field_timestep(*source, "temp")

fields = {}
for name in ["temp", "pres", "qv", "qi", "qc", "qs", "qr", "qg", "z_ifc", "u", "v", "w"]:
    nlevel = yac.get_field_collection_size(*source, name)
    fields[name] = Field.create(name, comp, points, nlevel,
                                dt, TimeUnit.ISO_FORMAT)
    yac.def_couple(*source, name,
                   "pamtra", "pamtra-grid", name,
                   dt, TimeUnit.ISO_FORMAT, Reduction.TIME_NONE,
                   interp)

yac.enddef()

with open(args.descriptorfile, "r") as dfile:
    exec(dfile.read())  # adds a variable "descriptorFile"


def call_pamtra(temp, pres, qv, qi, qc, qs, qr, qg, z_ifc, u, v, w,
                datetime):
    pamtra = pyPamtra.pyPamtra()
    for hydro in descriptorFile:
        pamtra.df.addHydrometeor(hydro)
    pamtra.nmlSet["passive"] = False
    pamtra.nmlSet["active"] = True
    pamtra.nmlSet["radar_mode"] = "simple"
    rh = 100*pyPamtra.meteoSI.q2rh(qv, temp, pres)
    pamtra.createProfile(hgt_lev=z_ifc.T[None, :, ::-1],
                         press=pres.T[None, :, ::-1],
                         temp=temp.T[None, :, ::-1],
                         relhum=rh.T[None, :, ::-1],
                         hydro_q=np.stack([qc.T, qi.T, qr.T, qs.T, qg.T], axis=-1)[:, ::-1, :],
                         )
    #pamtra.runPamtra(args.frequency)

    return pamtra # return the entire pamtra object and let the caller decide what to take as results

# prepare output
Ntimesperoutput = 360 #360 # every hour # TODO should be user config...

pam_time_idx = xr.IndexVariable(
        dims='time_idx',
        data=np.arange(0, Ntimesperoutput),
        attrs={'name':'time_idx',
               'long_name':'index of the time step within the output file'})

pam_hgt_idx = xr.IndexVariable(
        dims='hgt_idx',
        data=np.arange(0, nlevel-1),
        attrs={'name':'hgt_idx',
               'long_name':'height indx above surface',
               'direction':'bottom-up',
               'comment':'look at variable height to get the height profile above each location'})

pam_lev_idx = xr.IndexVariable(
        dims='lev_idx',
        data=np.arange(0, nlevel),
        attrs={'name':'lev_idx',
               'long_name':'level indx above surface',
               'direction':'bottom-up',
               'comment':'look at variable height to get the height profile above each location'})

pam_time = xr.DataArray(
        dims=['time_idx'],
        coords={'time_idx':pam_time_idx},
        data=np.arange(0, Ntimesperoutput),
        attrs={'name':'time',
               'long_name':'unixtime nanoseconds since 1970',
               'units':'nanoseconds'})

pam_lat = xr.DataArray(
        dims=['time_idx'],
        coords={'time_idx':pam_time_idx},
        data=np.empty((Ntimesperoutput,)),
        attrs={'name':'lat',
               'long_name':'location latitude',
               'unit':'degrees North'})

pam_lon = xr.DataArray(
        dims=['time_idx'],
        coords={'time_idx':pam_time_idx},
        data=np.empty((Ntimesperoutput,)),
        attrs={'name':'lon',
               'long_name':'location longitude',
               'unit':'degrees East'})

pam_hgt = xr.DataArray(
        dims=['hgt_idx', 'time_idx'],
        coords={'hgt_idx':pam_hgt_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel-1, Ntimesperoutput)),
        attrs={'name':'hgt',
               'long_name':'height above surface',
               'units':'meters'})

pam_Ze = xr.DataArray(
        dims=['hgt_idx', 'time_idx'],
        coords={'hgt_idx':pam_hgt_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel-1, Ntimesperoutput)),
        attrs={'name':'Ze',
               'long_name':'unattenuated reflectivity',
               'units':'dBZ'})

#################################
#"temp", "pres", "qv", "qi", "qc", "qs", "qr", "qg", "z_ifc", "u", "v", "w"
pam_w = xr.DataArray(
        dims=['lev_idx', 'time_idx'],
        coords={'lev_idx':pam_lev_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel, Ntimesperoutput)),
        attrs={'name':'w',
               'long_name':'wind vertical velocity',
               'units':'m/s'})
pam_v = xr.DataArray(
        dims=['hgt_idx', 'time_idx'],
        coords={'hgt_idx':pam_hgt_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel-1, Ntimesperoutput)),
        attrs={'name':'v',
               'long_name':'wind meridional velocity',
               'units':'m/s'})
pam_u = xr.DataArray(
        dims=['hgt_idx', 'time_idx'],
        coords={'hgt_idx':pam_hgt_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel-1, Ntimesperoutput)),
        attrs={'name':'u',
               'long_name':'wind zonal velocity',
               'units':'m/s'})
pam_z_ifc = xr.DataArray(
        dims=['lev_idx', 'time_idx'],
        coords={'lev_idx':pam_lev_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel, Ntimesperoutput)),
        attrs={'name':'z_ifc',
               'long_name':'level height above sea level',
               'units':'m'})
pam_qg = xr.DataArray(
        dims=['hgt_idx', 'time_idx'],
        coords={'hgt_idx':pam_hgt_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel-1, Ntimesperoutput)),
        attrs={'name':'qg',
               'long_name':'graupel mixing ratio',
               'units':'kg/kg'})
pam_qr = xr.DataArray(
        dims=['hgt_idx', 'time_idx'],
        coords={'hgt_idx':pam_hgt_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel-1, Ntimesperoutput)),
        attrs={'name':'qr',
               'long_name':'rain mixing ratio',
               'units':'kg/kg'})


pam_qs = xr.DataArray(
        dims=['hgt_idx', 'time_idx'],
        coords={'hgt_idx':pam_hgt_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel-1, Ntimesperoutput)),
        attrs={'name':'qs',
               'long_name':'snow mixing ratio',
               'units':'kg/kg'})

pam_qc = xr.DataArray(
        dims=['hgt_idx', 'time_idx'],
        coords={'hgt_idx':pam_hgt_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel-1, Ntimesperoutput)),
        attrs={'name':'qc',
               'long_name':'cloud drops mixing ratio',
               'units':'kg/kg'})

pam_qi = xr.DataArray(
        dims=['hgt_idx', 'time_idx'],
        coords={'hgt_idx':pam_hgt_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel-1, Ntimesperoutput)),
        attrs={'name':'qc',
               'long_name':'cloud ice mixing ratio',
               'units':'kg/kg'})

pam_qv = xr.DataArray(
        dims=['hgt_idx', 'time_idx'],
        coords={'hgt_idx':pam_hgt_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel-1, Ntimesperoutput)),
        attrs={'name':'qv',
               'long_name':'water vapour mixing ratio',
               'units':'kg/kg'})

pam_pres = xr.DataArray(
        dims=['hgt_idx', 'time_idx'],
        coords={'hgt_idx':pam_hgt_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel-1, Ntimesperoutput)),
        attrs={'name':'pres',
               'long_name':'air pressure',
               'units':'Pa'})

pam_temp = xr.DataArray(
        dims=['hgt_idx', 'time_idx'],
        coords={'hgt_idx':pam_hgt_idx,
                'time_idx':pam_time_idx},
        data=np.empty((nlevel-1, Ntimesperoutput)),
        attrs={'name':'temp',
            'long_name':'air temperature',
               'units':'K'})

#################################
variables = {'height':pam_hgt,
             'time':pam_time,
             'lat':pam_lat,
             'lon':pam_lon,
             'Ze':pam_Ze,
             'w':pam_w,
             'v':pam_v,
             'u':pam_u,
             'z_ifc':pam_z_ifc,
             'qg':pam_qg,
             'qr':pam_qr,
             'qs':pam_qs,
             'qc':pam_qc,
             'qi':pam_qi,
             'qv':pam_qv,
             'pres':pam_pres,
             'temp':pam_temp}

coordinates = {'time_idx':pam_time_idx,
               'hgt_idx':pam_hgt_idx}

import os
import socket
from datetime import datetime as datetimelib ### ugly but datetime is reused later and can be uglier
global_attributes = {'created_by':os.environ['USER'],
                     'host_machine':socket.gethostname(),
                     'created_on':str(datetimelib.now()),
                     'comment':'This is a test output of a aircraft'}

dataset = xr.Dataset(data_vars=variables,
                     coords=coordinates,
                     attrs=global_attributes)

itime = 0
while True:
    datetime = fields["temp"].datetime
    data = {}

    # TODO if datetime < first time or  > last time do nothing?
    timestep = pd.to_datetime(datetime)
    point_idx = list(track.time.values).index(track.sel(time=timestep, method='nearest').time) # TODO if point outside domain do nothing
    #print("point index {} {}".format(point_idx, type(point_idx)), file=sys.stderr)
    #print("lat {}, lon {} ".format(lats[point_idx], lons[point_idx]), file=sys.stderr)
    pam_time[itime] = timestep

    for v, f in fields.items():
        varfield, info = f.get()
        datav = varfield[..., point_idx] # we want this to be (Nhgt, something)
        if datav.ndim == 1:
            data[v] = datav[:, np.newaxis]
        else:
            data[v] = datav

    pam_lat[itime] = np.rad2deg(lats[point_idx])
    pam_lon[itime] = np.rad2deg(lons[point_idx])


    pam_z_ifc[:, itime] = data['z_ifc'].T[0, ::-1]
    pam_w[:, itime] = data['w'].T[0, ::-1]
    pam_u[:, itime] = data['u'].T[0, ::-1]
    pam_v[:, itime] = data['v'].T[0, ::-1]
    pam_qg[:, itime] = data['qg'].T[0, ::-1]
    pam_qr[:, itime] = data['qr'].T[0, ::-1]
    pam_qs[:, itime] = data['qs'].T[0, ::-1]
    pam_qc[:, itime] = data['qc'].T[0, ::-1]
    pam_qi[:, itime] = data['qi'].T[0, ::-1]
    pam_qv[:, itime] = data['qv'].T[0, ::-1]

    pam_temp[:, itime] = data['temp'].T[0, ::-1]
    pam_pres[:, itime] = data['pres'].T[0, ::-1]

    #pamtra = call_pamtra(**data, datetime=datetime) # this function does not actually run pamtra
    #pam_Ze[:, itime] = pamtra.r['Ze'][0,0,:,0,0,0] # this extracts only 1 frequency... need to think about it
    #pam_hgt[:, itime] = pamtra.r['radar_hgt'][0, 0, :]


    itime += 1
    if itime == Ntimesperoutput:
        dataset.to_netcdf('aircraft_{}.nc'.format(pd.to_datetime(pam_time[0].values).strftime('%Y%m%d-%H%M%S')))
        itime = 0

    if info == Action.GET_FOR_RESTART:
        break

