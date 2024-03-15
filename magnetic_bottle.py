
from math import *
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

####################################################

class problem():
    def __init__(self): ## Initialize instance with the charge of each particle ## charges are pre-set for now for simplicity
        self.particle_charge = 1
        self.particle_mass = 1

    def v_prim(self, direction, vector_v, vector_b): ## check whether input data 
        if(direction == 'x'):
            return (self.particle_charge/self.particle_mass)*(vector_v['y']*vector_b['z'] - vector_v['z']*vector_b['y']) ## solving for acceleration of particle in x vector
        if(direction == 'y'):
            return (self.particle_charge/self.particle_mass)*(vector_v['z']*vector_b['x'] - vector_v['x']*vector_b['z']) ## solving for acceleration of particle in y vector
        if(direction == 'z'):
            return (self.particle_charge/self.particle_mass)*(vector_v['x']*vector_b['y'] - vector_v['y']*vector_b['x']) ## solving for acceleration of particle in z vector
    def particle_pos(self, current_pos, v_prim, step_size):

        for component in current_pos:
            current_pos[component] = current_pos[component] + (v_prim[component] * step_size)

        return current_pos 

class solver():

    def __init__(self, step_size, target_x, problem): ## artifically made data used to start off the particles motion

        self.vector_pos = {'x' : 1, 'y' : -1, 'z' : 0.001} ### position of particle on catersian plane (x,y,z)
        self.vector_v = {'x' : (1-(0.5*0.5)), 'y' : (1-(0.5*0.5)), 'z' : 0.2} ## 
        self.vector_v_next = {'x' : None, 'y' : None, 'z' : None}
        self.vector_b = {'x' : 0, 'y' : 0, 'z' : 15}

        self.step_size = step_size
        self.target_x = target_x
        self.handles = {}
        self.handles_read = {}
        self.problem = problem
        self.debug = False

######## File Saving fn Start Here ########
        

    def load_handle_write(self,file_name):
        if not file_name in self.handles:
            if(os.path.exists(file_name)):
                os.remove(file_name)
            #Open file/files for saving
            self.handles[file_name] = open(file_name, 'wb')
        return self.handles[file_name]

    def load_handle_read(self,file_name):
        if not file_name in self.handles_read:
            #Open file/files for saving
            self.handles_read[file_name] = open(file_name, 'rb')
        return self.handles_read[file_name]

    def load(self, file_name, debug):
        if(debug and self.debug):
            print("Loading values from:", file_name)
        values = self.load_handle_read(file_name + '.txt').readlines()
        self.load_handle_read(file_name + '.txt').seek(0)
        for i, value in enumerate(values):
            # Decode each bytes object to a string before stripping and converting
            values[i] = float(value.decode('utf-8').rstrip('\n'))
        return values
    
    def save(self, file_name, value, debug):
    
        if(debug and self.debug):
            print("Saving to file [\"" + debug + "\"]:", value)
        self.load_handle_write(file_name + '.txt').write((str(value) + "\n").encode('utf-8'))


######## End of File Saving fn ########
        
    def euler(self):
        
        #No of steps to take!
        no_of_steps = int(self.target_x/self.step_size + 2) # +1 for range, +1 for getting the division right :)
        
        print("\nRunning " + str(no_of_steps) + " no. of steps!")
        
        # Set m
        m = 15
        c_value = 9

        for i in range(no_of_steps):
            
            if i < 1:
            
                #First step doesn’t need to be calculated
                pass
            
            if i > 0:
            
                # Update magnetic field lines to create "convergence" effect. This should eventually let the particle "flip" direction as well.
                # Bz = -((z - c)^2) + c^2
                self.vector_b['z'] = -((self.vector_pos['z'] - c_value)*(self.vector_pos['z'] - c_value)) + c_value*c_value ## line 107
                if self.vector_b['z']<0: 
                    self.vector_b['z'] = 0
                    self.vector_b['x'] = 0
                    self.vector_b['y'] = 0
                else:
                    # Bx = -(1/2) * m * x * (((z - c)^2) / c^2) / (x^2 + y^2)
                    self.vector_b['x'] = -(0.5) * (m * self.vector_pos['x'])*(((self.vector_pos['z'] - c_value)*(self.vector_pos['z'] - c_value))/(c_value*c_value))/(self.vector_pos['x']*self.vector_pos['x'] + self.vector_pos['y']*self.vector_pos['y'])
                    # By = -(1/2) * m * y * (((z - c)^2) / c^2) / (x^2 + y^2)
                    self.vector_b['y'] = -(0.5) * (m * self.vector_pos['y'])*(((self.vector_pos['z'] - c_value)*(self.vector_pos['z'] - c_value))/(c_value*c_value))/(self.vector_pos['x']*self.vector_pos['x'] + self.vector_pos['y']*self.vector_pos['y'])

                # For each component in the speed, update it using acceleration! :)
                for component in self.vector_v:
                
                    # Calculate v value from it’s derivative with good old euler :)
                    acceleration = self.problem.v_prim(component, self.vector_v, self.vector_b) * self.step_size
                    self.vector_v_next[component] = self.vector_v[component] + acceleration

                # Update particle position based on velocity and previous position.
                self.vector_pos = self.problem.particle_pos(self.vector_pos, self.vector_v, self.step_size)

                # Save all known data for future plotting, this could be rewritten with a "listener"/"server" system to minimize time
                # spent writing data to disk.
                self.save('part_speed_x', self.vector_v['x'], 'speed_x')
                self.save('part_speed_y', self.vector_v['y'], 'speed_y')
                self.save('part_speed_z', self.vector_v['z'], 'speed_z')
                                                        
                self.save('part_pos_x', self.vector_pos['x'], 'pos_x')
                self.save('part_pos_y', self.vector_pos['y'], 'pos_y')
                self.save('part_pos_z', self.vector_pos['z'], 'pos_z')
                self.save('time', i, 'time')
                self.save('magnetic_vector_bz', self.vector_b['z'], 'magnetic_z')
                self.save('magnetic_vector_bx', self.vector_b['x'], 'magnetic_x')
                self.save('magnetic_vector_by', self.vector_b['y'], 'magnetic_y')
                
                # Update old speed to new speed.
                self.vector_v = self.vector_v_next
           
        for i, handle in self.handles.items():
            handle.close()

######## Plotting fn Start Here ########
            

    def plot_show(self):
        #Let’s show the graph as well
        plt.show()

    def paramteric_curve(*args):
        #Initialize a new graph
        plt.figure()
        
        #Need at least x and y values to begin!
        if(len(args) < 3):
            print("Error - please provide at least two variables to plot!")
            sys.exit()

        #Manually define self since we are running with a variable no of variables!
        self = args[0]
        x_val = self.load(args[1], True)
        
        for i, arg in enumerate(args):
            if i > 1:
                values = self.load(arg, True)

                #Debug info
                if(self.debug):
                    for i, x in enumerate(x_val):
                        print("Plotting: ", x_val[i], values[i])
                
                #Plot each value-set
                plt.plot(x_val, values, '.-')
                plt.xlabel(args[1])
                plt.ylabel(arg)
                plt.title(arg + " vs " + args[1])
                plt.grid(True)
        
    #def plot_3d(*args):

    def plot(*args):
        #Initialize a new graph
        plt.figure()
        
        #Need at least x and y values to begin!
        if(len(args) < 3):
            print("Error - please provide at least two variables to plot!")
            sys.exit()

        #Manually define self since we are running with a variable no of variables!
        self = args[0]
        x_val = self.load(args[1], True)
        
        for i, arg in enumerate(args):
            if i > 1:
                values = self.load(arg, True)

                #Debug info
                if(self.debug):
                    for i, x in enumerate(x_val):
                        print("Plotting: ", x_val[i], values[i])
                
                #Plot each value-set
                plt.plot(x_val, values, '.-')
                plt.xlabel(args[1])
                plt.ylabel(arg)
                plt.title(arg + " vs " + args[1])
                plt.grid(True)
                

######## End of Plotting fn ######## 
    def custom_vars(*args):
        self = args[0]

        if(len(args) < 4):
            sys.exit('Error')
        
        # We are given variable names, load them from the files
        pass_args = []

        for arg in args[3:]:
            pass_args.append(self.load(arg, True))
            
        for var_no, var in enumerate(self.load(args[3], True)):
            
            # We only want to give one value per "variable" on in each iteration
            temp_pass_args = []
        
            for arg_no, arg in enumerate(pass_args):
            
                temp_pass_args.append(pass_args[arg_no][var_no])
            
            # Run the custom function
            self.save(args[1], args[2](*temp_pass_args), 'Adding custom value')
        
        for i, handle in self.handles.items():
            handle.close()

#Solve this using two euler approximations at the same time. Predict first derivative then second using the first.
#Initial values - x = 0, x’ = 1, x’’ = 1
# Usage of custom functions:
# 1) Define custom function with variables you want from previous data
# 2) Run custom_vars(), with arg1 -> name of variable you want to save it to, arg2 -> name of custom function,
# arg3 -> the "x" value you want to iterate with, can be any value, arg(3+n) -> any other variables you
# want to use in your custom function

#Config
step_size = 0.001 #Accuracy - "stepsize" - adjust for more/less precise plots, experiment to trade-off time/cpu-hours#Run it

problem = problem() #Define our problem

#Run it for several stepsizes

print("\nDoing stepsize:", float(0.2))
solver_instance = solver(step_size, 20, problem) #Config our solver
solver_instance.euler()
solver_instance.plot('time', 'part_speed_z', 'magnetic_vector_by', 'magnetic_vector_bx') #Plot the speed versus the z pos
solver_instance.plot('part_pos_z', 'magnetic_vector_bz') #Plot the speed versus the z pos

#Close all open handles after we’ve shown the graphs!

print("Here is vector_b:", '\n', solver_instance.vector_b)
print("Here is vector_v:", '\n', solver_instance.vector_v)

try:
    solver_instance.plot_show()
except:
    print("Quitting!")
    for i, handle in solver_instance.handles.items():
        handle.close()
    for i, handle in solver_instance.handles_read.items():

        handle.close()