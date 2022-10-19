# fmesh

Simple mesh generation for regional FESOM2 setups.

Originally created by Sergey Kirillov.

Video with process of defining regions and looking at the result can be found under the link below: 
[<img src="https://user-images.githubusercontent.com/3407313/196788021-4a86fc91-0883-4e8d-8e4c-21c0b6bcdb65.png" width="50%">](https://youtu.be/8sxx6JwoLsg "Creating regional FESOM2 mesh with fmesh.")

## Installation

Create new `conda` environment from `requirements.yml` by executing:

```bash
conda env create -f requirements.yml
```
and activate it
```
conda activate fmesh
```

#### Downloading necessary data

You have to download bathymetry and coastlines. This can be done by executing:

```bash
./download_data.sh
```

## Usage

You have to configure your future mesh in the `configure.yaml` and then run:

```bash
python fmesh.py
```
in the `fmesh` directory.

## Configuration

The process of mesh creation is split into several steps, that you can configure. In general you define the base resolution (ether constant, or increasing towards poles) and perform refinement on some regions of the basis mesh.

### Define base resolution

The main step is defining your base resolution. There are two options here:
- constant resolution (`do_mercator_refinement: false`)
- resolution changing with latitude (`do_mercator_refinement: true`)

For the resolution changing with latitude you can define the latitude where resolution become constant ( `norhtern_boundary`, `southern_boundary`) and the value of constant resolution (`norhtern_lowest_resolution`, `southern_lowest_resolution`).

### Refine resolution along coastlines

To turn it on you have to set `do_refine_along_coastlines: true`. You can control:
- `min_resolution` resolution at the coast.
- `max_distance` distance from the coast, where resolution is set back to base resolution. So the resolution will change from `min_resolution` at the coast to base mesh resolution at `max_distance`.
- `min_length` minimum length of the coastline to be considered. For example if it set to 200, the island with 200 km coastline will not be considered.
- `averaging` the smothing lenghth of the coastline.

### Select regions with increased resolution

Regions are defined in .kml files, that contain two paths for inner region (where the resolution will be set to constant `resolution` value) and outer region, where resolution will be set to the base resolution. In the area between the paths resolution will be linearly changing from `resolution` value to base resolution.

The easiest way to create paths is to use Goole Earth, and path tool. The paths should not be closed, this will be done automatically by `fmesh`.
<img width="1607" alt="Screenshot 2022-10-19 at 21 29 04" src="https://user-images.githubusercontent.com/3407313/196785898-251c66ee-5822-4042-b157-f66318d09962.png">

Additional parameters for each region:
- `precision` each section of polygon will by splitted in `precision` nuber of times, to make linear interpolation more precise.
- `order` If True, fist polygon is internal and second is external. I False, then other way round.

### Set vertical levels

For convinience you can set vertical levels in your resulting `aux3d.out` file with `levels` parameter.

### First look at the mesh

The `fmesh` also saves .vtk file, with resulting mesh, that can be opened in [ParaView](https://www.paraview.org/):
<img width="1920" alt="Screenshot 2022-10-19 at 21 40 35" src="https://user-images.githubusercontent.com/3407313/196788021-4a86fc91-0883-4e8d-8e4c-21c0b6bcdb65.png">

 
