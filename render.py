from direct.showbase.ShowBase import ShowBase
from math import isnan
from os.path import abspath
from panda3d.core import Filename
from particles import c, gravity, m, MAX, panda, r, v
from scipy.constants import pi, Planck
from scipy.special import betainc, gamma
from sys import path
clear, con = [0] * len(c), [[False] * len(c)] * len(c)
def hist(i, f): pass # Add i to the frame before f iterations
def rerender(f): pass # Rerender the last f frames because of time travel and create a parallel universe
def norm(l, r): return sum([(l[i] - r[i]) ** 2 for i in range(len(l))]) ** 0.5
for i in range(len(c)):
    for j in range(len(c)): con[i][j] = norm(c[i], c[j]) < r[i] + r[j]
class View(ShowBase):
    def __init__(self):
        ShowBase.__init__(self)
        self.taskMgr.add(self.update, "UpdateScene")
        self.balls, self.dir = [], Filename.fromOsSpecific(abspath(path[0])).getFullpath()
        for i in range(len(c)):
            if panda: self.balls += [Actor('/models/panda-model')]
            else:
                self.balls += [loader.loadModel(self.dir + '/sphere.egg')]
                if not m[i]: self.balls[-1].setColor(255, 255, 0, 1) # Show photons in yellow
                if m[i] > 0: self.balls[-1].setColor(255, 255, 255, 1) # Show masses in white
                if m[i] < 0: self.balls[-1].setColor(0, 0, 0, 1) # Show negative masses in black
            self.balls[-1].setScale(r[i], r[i], r[i])
            self.balls[-1].setPos(c[i][0], c[i][1], c[i][2])
            self.balls[-1].reparentTo(self.render)
    def update(self, task):
        global c, m, v, r, con, clear
        TimeDistortion = 0
        for i in range(len(c)):
            if m[i] > 0: #apply the lorentz transform, which is not applied when kasimir effect takes place
                for p in range(len(v[i])): v[i][p] *= (1 - (v[i][p] / (3 * 10 ** 8)) * 2) ** 0.5
            normC = c[i] if not(any(c[i])) else [p / (norm(c[i], [0] * len(c[i]))) for p in c[i]]
            try:
                clear[i] = int((1 - (norm(v[i], [0] * len(v[i])) / (3 * 10 ** 8)) ** 2) **- 0.5) #time dialation
                self.balls[i].hide()
            except TypeError: pass # that means the time dilation is complex so this is the kasimir effect which will be fixed in the hist function
            if not(any(v[i])) or Planck / norm(v[i], [0] * len(v[i])) / m[i]: #de broglie wave length is smaller than 6.2*10^34Hz therefore the wave length is smaller than plank's constant
                clear[i] = 2
                self.balls[i].hide()
            if m[i] <= 0 and norm(v[i], [0] * len(v[i])) > 3 * 10 ** 8: #if a particle passed the speed of light a reverse time dialation is applied
                d = int((1 - ((3 * 10 ** 8) / norm(v[i], [0] * len(v[i]))) ** 2) **- 0.5)
                hist(i, d)
                TimeDistortion = max(TimeDistortion, d)
            for j in range(len(c)):
                if i == j: continue
                dis = norm(c[i], c[j])
                if dis > r[i] + r[j]: #if there is no contact the effective force is gravity
                    con[i][j], con[i][j] = False, False
                    for p in range(len(c[i])): v[j][p] += 666 * normC[p] * m[i] / (dis ** 2) / (10 ** 13)
                elif not con[i][j]: #if this is the initial contact the effective force is elastic collision
                    con[i][j], con[i][j] = True, True
                    for p in range(c[i]):
                        vj, vi = v[j][p], v[i][p]
                        v[j][p], v[i][p] = 2 * m[i] * vi / (m[i] + m[j]) + (m[j] - m[i]) * vj / (m[i] + m[j]), 2 * m[j] * vj / (m[j] + m[i]) + (m[i] - m[j]) * vi / (m[j] + m[i])
                else: #if there is already initial contact the effective force is multidimentional drag
                    for p in range(c[i]):
                        h = r[i] + r[j] - dis
                        if r[j]: v[j][p] += gamma(len(c[i]) + 1) * r[j] ** (len(c[i]) - 1) * betainc(4 * h ** 2 - 4 * h ** 3 / r[j] + h ** 4 / r[j] ** 2, (len(c[i]) - 1) / 2, 0.5) * m[i] * (v[i][p] - v[j][p]) ** 2 / (4 * gamma(len(c[i]) / 2) * r[i] ** len(c[i]))
        for i in range(len(c)):
            if clear[i]:
                clear[i] -= 1
                if clear[i] == 1: self.balls[i].show()
            if gravity: v[i][1] -= 9.8
            for p in range(len(c[i])):
                c[i][p] = max(min(c[i][p] + v[i][p].real, MAX - c[i][p] + v[i][p].real), c[i][p] + v[i][p].real - MAX)
                if isnan(c[i][p]): raise TypeError("Need to be fixed")
        if TimeDistortion: rerender(TimeDistortion) #if the space-time has changed do stuff
        for i in range(len(c)): self.balls[i].setPos(c[i][0], c[i][1], c[i][2])
        return task.cont
View().run()
#TODO: write hist and rerender