# Background and Validation

This wiki is built in Notion. Here are all the tips you need to contribute.

# Installation

To install using pip Run:

```python
python3 -m pip install 2D_Panel CFD
```

Or:

Download source folder anywhere in your machine & set it to PYTHONPATH

```python
export PYTHONPATH="${PYTHONPATH}:/path/to/source_directory/src"
```

It's advisable to run this package from [RUN-spyder.py](http://RUN-spyder.py) through an IDE like spyder for ease of use, and prolonged variable storage.

## Executables

<aside>
💡 The source directory should be set up as PYTHONPATH

</aside>

- **RUN.py** is meant to be run from terminal.
- **RUN_spyder.py** can be run with an IDE, such as spyder to improve multi-Core Utilisation & post-processing experience. ****
- The prompts makes the program easy to run, & visualise for any 2D sample problem.

# General Background

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/README/Map-1_Step-.jpg" alt="Flow over a cylinder" width="1500"/>

Flow over a cylinder

---

> **The project has been started as a Open Source repository for CFD solvers.**
> 

> **This version comprises of a 2D Staggered Grid with Inlet, Outlet & Wall Boundary conditions. Obstacles can be imported & transformed with a list of points or with the inbuilt elliptical geometries.**
> 

> **First order Upwind Scheme is used for Velocity with very good results for the benchmark Lid Driven Cavity problem when compared to results in Ghia etal.**
> 

> **The SIM runs stable with Python for <10000 Cells after which Residual plotting becomes laggy. spyder (Anaconda IDE) provides great speed-ups with multi-core utilisation & also improves the post-processing experience. Some minor modifications in concurrent Residual plotting makes running SIM in spyder a better solution for now. The Sequential prompts based model is based on a GUI approach and will be ported to it in the next update.**
> 

> **The lack of multi-threading support in python trumps the ease of accessibility of matplotlib library. We will be looking to port into C++ immediately utilizing vtk libraries with paraview & blender for visualization.**
> 

> **The framework is designed to test new FVM schemes, & Coupling solvers. All popular convection schemes will be added soon. Multiple solvers will be available in the next updates, the likes of SIMPLER, PISO, Pimple etc. Future plans also include Unsteady & VOF solvers.**
> 

> **The program works as a sequential prompt, for SIM Parameters. The prompts are designed keeping in mind a GUI approach, which will be available in the next update. There are frequent Check Cycles to render the result & modify any inputs. We'll go through an exemplary First Run in the next Section.**
> 

# Validation of Solver

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/README/ezgif.com-gif-maker(3).gif" alt="Vortex Shedding flow over a cylinder" width="1500"/>


Vortex Shedding flow over a cylinder

---

> For validation of the solver laid out, following strategies are used:
> 
1. Comparison with Benchmark Problem Lid Driven Cavity 
    1. Reference study Ghia etal. Re = 100, 1000, 5000

## Lid Driven Cavity Benchmark Ghia etal.

### **Residuals**

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/README/Untitled.png" alt="Untitled" width="1500"/>

### **Benchmark Test at Re=100**

- First Order Upwind scheme

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/README/Untitled%201.png" alt="Untitled" width="1500"/>

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/README/Untitled%202.png" alt="Untitled" width="1500"/>

### **Benchmark Test at Re=1000**

- First Order Upwind scheme

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/README/Untitled%203.png" alt="Untitled" width="1500"/>

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/README/Untitled%204.png" alt="Untitled" width="1500"/>

### **Benchmark Test at Re=5000**

- First Order Upwind scheme

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/README/Untitled%205.png" alt="Untitled" width="1500"/>

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/README/Untitled%206.png" alt="Untitled" width="1500"/>


### Conclusion

First order UPWIND Scheme is good for low Reynolds no. but is only first order accurate to capture higher gradient. 

## Fully developed flow between Parallel Plates

### Velocity Profile [at X=Lx/2 and Y=0.8Ly]

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/README/Untitled%207.png" alt="Untitled" width="1500"/>

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/README/Untitled%208.png" alt="Untitled" width="1500"/>

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/README/Map-1_Step-200.jpg" alt="Map-1 Step-[200].jpg" width="1500"/>

### Conclusion

First order UPWIND Scheme with high y-gradient.
