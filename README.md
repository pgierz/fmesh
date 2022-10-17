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
- constant resolution 
