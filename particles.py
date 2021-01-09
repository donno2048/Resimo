MAX = 100 # The radius of the "world"
c = [[0, 0, 0], [1, 2, 3]] # Initial coordinates (it works even in higher dimentions...)
m = [1, 4] # Masses of the balls (can be negative thanks to the casimir effect)
v = [[0, 0, 0], [0, 0, 1]] # Initial velocity coordinates
r = [2, 1] # Balls' radius
gravity = False # Wheter to apply gravity
panda = False # Use pandas instead of spheres ;)