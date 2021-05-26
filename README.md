# Automatic regulation for OpenFOAM

This is part of a masters thesis: "Automatic Regulation of Industrial Processes in OpenFOAM"

This project aims to provide a way to simulate automatic regulation in OpenFOAM CFD framework. It's two main components are:
- `Regulator`, a component responsible for generating control signal based on current system state and specified target value.
- Boundary conditions using `Regulator`, acting as an acutor and used to integrate the Regulator with OpenFOAM case.

## Examples
### 1. Laminar flow in pipe
Task: control flow rate at inlet to maintain target temperature at outlet. PID algotirhm is used.

### 3. Temperature regulation in room
Task: control heat flux from floor to maintain a target temperature of the room. On/off regulation is used.
