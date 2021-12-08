# Background and Validation

This wiki is built in Notion. Here are all the tips you need to contribute.

# General Background

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/img/Map-1_Step-.jpg" alt="Flow over a cylinder" width="700"/>

Flow over a cylinder

---

> **The project has been started as a Open Source repository for CFD solvers. The motive is to provide handy easy to understand code with multitude of CFD schemes for cfd developers. Also, needs to remain functional as an easy to setup open source solver for users. This release only comprises of a terminal sequential prompt, simple and effective. We have immediate plans of implementing a PyQT GUI to it.**
> 

> Head to the notion page for more information on how to add to this project:
> 

[https://florentine-hero-1e6.notion.site/2D_Panel-CFD-ad63baa924ee4a32af8a52b8134c0360](https://www.notion.so/2D_Panel-CFD-ad63baa924ee4a32af8a52b8134c0360)

> **This version comprises of a 2D Staggered Grid with Inlet, Outlet & Wall Boundary conditions. Obstacles can be imported & transformed with a list of points or with the inbuilt elliptical geometries.**
> 

> **First order Upwind Scheme is used for Velocity with very good results for the benchmark Lid Driven Cavity problem when compared to results in Ghia etal.**
> 

> **The SIM runs stable with terminal-python for <10000 Cells after which Residual plotting becomes laggy. spyder (Anaconda IDE) provides great speed-ups with multi-core utilisation & also improves the post-processing experience. The Sequential prompts based model is based on a GUI approach and will be ported to it in the next update.**
> 

> **The lack of multi-threading support in python trumps the ease of accessibility of matplotlib library. We will be looking to port into C++ immediately utilizing vtk libraries with paraview & blender for visualization.**
> 

> **The framework is designed to test new FVM schemes, & Coupling solvers. All popular convection schemes will be added soon. Multiple solvers will be available in the next updates, the likes of SIMPLER, PISO, Pimple etc. Future plans also include Unsteady & VOF solvers.**
> 

> **The program works as a sequential prompt, for SIM Parameters. The prompts are designed keeping in mind a GUI approach, which will be available in the next update. There are frequent Check Cycles to render the result & modify any inputs. We'll go through an exemplary First Run in the next Section.**
> 

# Installation

### Method: 1

To install using pip Run:

```python
python3 -m pip install 2D_Panel-CFD
```

Or:

### Method: 2

[https://github.com/Fluidentity/2D_Panel-CFD](https://github.com/Fluidentity/2D_Panel-CFD)

- Clone github `[RUN_package](https://github.com/Fluidentity/2D_Panel-CFD.git)` to anywhere in your machine from:

```tsx
cd /insert/folder/address/cfd
git clone https://github.com/Fluidentity/2D_Panel-CFD.git
```

- Set it to PYTHONPATH with:

```python
export PYTHONPATH="${PYTHONPATH}:/insert/folder/address/cfd/RUN_package"
```

It's advisable to run this package from [RUN-spyder.py](http://RUN-spyder.py) through an IDE like spyder for ease of use, and prolonged variable storage. Also, spyder has some great plotting interface.

## Executable

<aside>
💡 The source directory should be set up as PYTHONPATH if not installed using pip

</aside>

### Method: 1

Open python environment with: (in terminal)

```python
python3
```

or (if python —version is >3)

```python
python
```

then insert:

```python
from RUN_package import RUN
```

- **RUN.py** is meant to be run from terminal.

### Method: 2

Run on IDE by cloning [RUN_package](https://github.com/Fluidentity/2D_Panel-CFD/tree/main/RUN_package) from Github.

> Open python IDE like spyder from RUN_package directory:
> 

> Run RUN-spyder.py
> 

The cells for pre-processor, solver & post processors are different. Need to run all. 

- **RUN_spyder.py** can be run with an IDE, such as spyder to improve multi-Core Utilisation & post-processing experience. ****

# Validation of Solver

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/img/ezgif.com-gif-maker(3).gif" alt="Vortex Shedding flow over a cylinder" width="700"/>

Vortex Shedding flow over a cylinder

---

> For validation of the solver laid out, following strategies are used:
> 
1. Comparison with Benchmark Problem Lid Driven Cavity 
    1. Reference study Ghia etal. Re = 100, 1000, 5000

## Lid Driven Cavity Benchmark Ghia etal.

### **Residuals**

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/img/Untitled.png" alt="Untitled" width="400"/>

### **Benchmark Test at Re=100**

- First Order Upwind scheme

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/img/Untitled%201.png" alt="Untitled" width="400"/>

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/img/Untitled%202.png" alt="Untitled" width="400"/>

### **Benchmark Test at Re=1000**

- First Order Upwind scheme

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/img/Untitled%203.png" alt="Untitled" width="400"/>

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/img/Untitled%204.png" alt="Untitled" width="400"/>

### **Benchmark Test at Re=5000**

- First Order Upwind scheme

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/img/Untitled%205.png" alt="Untitled" width="400"/>

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/img/Untitled%206.png" alt="Untitled" width="400"/>


### Conclusion

First order UPWIND Scheme is good for low Reynolds no. but is only first order accurate to capture higher gradient. 

## Fully developed flow between Parallel Plates

### Velocity Profile [at X=0.8*Lx and Y=0.5*Ly]

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/img/Untitled%207.png" alt="Untitled" width="400"/>

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/img/Untitled%208.png" alt="Untitled" width="400"/>

<img src="https://github.com/Fluidentity/2D_Panel-CFD-/blob/main/img/Map-1_Step-200.jpg" alt="Map-1 Step-[200].jpg" width="1500"/>

### Conclusion

The Umax Velocity comes close to 1.5 feactor for steady flow between parallel plates. First order UPWIND Scheme with high y-gradient.



