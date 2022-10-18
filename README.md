# fmesh

Simple mesh generation for regional FESOM2 setups.

## Installation

Create new `conda` environment from `requirements.yml` by executing:

```bash
conda env create -f requirements.yml
```
and activate it
```
conda activate fmesh
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
- `min_resolution` resolution at the coast
- `max_distance` distance from the coast, where resolution is set back to base resolution. So the resolution will change from `min_resolution` at the coast to base mesh resolution at `max_distance`.
- `min_length` minimum length of the coastline to be considered. For example if it set to 200, the island with 190 km coastline will not be considered.
- `averaging` the smothing lenghth of the coastline.

 
