# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as patches
from matplotlib.path import Path
plt.rcParams["figure.figsize"] = (8,8)

def plot_orbital_boxes (ax, x, y, boxes_number, electrons_number, box_side = 1, spacing_f = 5):
    Xi = x - boxes_number*box_side/2.
    Yi = y - box_side/2.
    def add_spin(Xi,Yi,box_side,direction):
        unit = box_side*0.8
        spacing = unit/float(spacing_f)
        hspacing = spacing/2.
        v_pad = box_side*0.1 + Yi
        h_pad = box_side/2. + Xi

        codes = [Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
                 ]
        if direction  == 'down':
            verts_down = [
                (h_pad + hspacing, v_pad + unit ), # left, bottom
                (h_pad + hspacing, v_pad ), # left, top
                (h_pad + hspacing + unit/5.,v_pad + unit*0.40), # right, top
                (h_pad + hspacing + unit/20., v_pad + unit*0.40), # right, bottom
                (h_pad + hspacing + unit/20., v_pad + unit),
                (h_pad + hspacing, v_pad + unit), # ignored
                ]
            spin_down = patches.PathPatch(Path(verts_down, codes),
                                          facecolor='k',
                                          lw=0.1,
                                          zorder=10)
            return spin_down

        if direction == 'up':

            verts_up = [
                (h_pad - hspacing, v_pad ), # left, bottom
                (h_pad - hspacing, v_pad + unit), # left, top
                (h_pad - hspacing - unit/5., v_pad +unit*0.60), # right, top
                (h_pad - hspacing - unit/20,v_pad + unit*0.60), # right, bottom
                (h_pad - hspacing - unit/20,v_pad),
                (h_pad - hspacing, v_pad), # ignored
                ]

            spin_up = patches.PathPatch(Path(verts_up, codes),
                                  facecolor='k',
                                  lw=0.1,zorder=10)
            return spin_up


    # plot the spins using Aufbau
    if electrons_number > 0:
        moduloelectrons = electrons_number%boxes_number
        if moduloelectrons > boxes_number:
            Warning ("electrons_number greater than boxes number")
        if moduloelectrons == 0:
            moduloelectrons = 1
        if electrons_number <= boxes_number:
            for j in range(electrons_number):
                ax.add_patch(add_spin(Xi+box_side*j,Yi,box_side,direction='up'))
        else:
            for e in range(moduloelectrons):
                ax.add_patch(add_spin(Xi+box_side*e,Yi,box_side,direction='down'))
            for j in range(boxes_number):
                ax.add_patch(add_spin(Xi+box_side*j,Yi,box_side,direction='up'))

class ED:
    def __init__(self,energy_gap=10, aspect='equal'):
        # plot parameters
        self.ratio = 1.5
        self.dimension = 'auto'
        self.space = 'auto'
        self.offset = 'auto'
        self.offset_ratio = 0.1
        self.color_bottom_text = 'black'
        self.aspect = aspect
        self.energy_gap = energy_gap
        # data
        self.pos_number = 0
        self.energies = []
        self.positions = []
        self.colors = []
        self.top_texts = []
        self.bottom_texts = []
        self.left_texts = []
        self.right_texts = []
        self.links = []
        self.arrows = []
        self.electons_boxes = []
        # matplotlib fiugre handlers
        self.fig = None
        self.ax = None

    def add_level(self, energy, bottom_text='', position=None, color='k', top_text='', right_text='', left_text=''):

        if position is None:
            position = self.pos_number + 1
            self.pos_number += 1
        elif position == 'last':
            position = self.pos_number
        if top_text == 'Energy':
            top_text = energy

        self.colors.append(color)
        self.energies.append(energy)
        self.positions.append(position)
        self.top_texts.append(top_text)
        self.bottom_texts.append(bottom_text)
        self.left_texts.append(left_text)
        self.right_texts.append(right_text)
        self.links.append([])
        self.arrows.append([])


        # print ('self.positions')
        # print (self.positions)
        # print ('==============')

    def add_link(self, start_level_id, end_level_id, color='k',ls='--',linewidth=0.7,coeficient=0,reverse=True):
        self.links[start_level_id].append((end_level_id, ls, linewidth, color,coeficient,reverse))

        # print ('self.add_link')
        # print (self.links)
        # print ('==============')


    def add_electronbox(self,level_id,boxes,electrons,side=0.8,spacing_f=5):
        self.__auto_adjust(self.energy_gap)
        side= side*self.energy_gap*0.1
        x = self.positions[level_id]*(self.dimension+self.space)+self.dimension*0.5
        y = self.energies[level_id]
        self.electons_boxes.append((x, y, boxes, electrons, side, spacing_f))


    # def plot(self, show_IDs=False,ylabel="Energy / $kcal$ $mol^{-1}$"):
    def plot(self, show_IDs=False,ylabel="Energy / $eV$"):
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect=self.aspect)
        ax.set_ylabel(ylabel)
        ax.axes.get_xaxis().set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        self.__auto_adjust(self.energy_gap)

        data = zip(self.energies,  # 0
                   self.positions,  # 1
                   self.bottom_texts,  # 2
                   self.top_texts,  # 3
                   self.colors,  # 4
                   self.right_texts, # 5
                   self.left_texts,)  # 6

        for level in data:
            start = level[1]*(self.dimension+self.space)
            ax.hlines(level[0], start, start + self.dimension, color=level[4])
            ax.text(start+self.dimension/2.,  # X
                    level[0]+self.offset*0.4,  # Y
                    level[3],  # self.top_texts
                    horizontalalignment='center',
                    verticalalignment='bottom')

            ax.text(start + 1.1*self.dimension,  # X
                    level[0],  # Y
                    level[5],  # self.right_text
                    horizontalalignment='left',
                    verticalalignment='center',
                    color=self.color_bottom_text)

            ax.text(start - 0.1*self.dimension,  # X
                    level[0],  # Y
                    level[6],  # self.bottom_text
                    horizontalalignment='right',
                    verticalalignment='center',
                    color=self.color_bottom_text)

            ax.text(start + self.dimension/2.,  # X
                    level[0] - self.offset*2,  # Y
                    level[2],  # self.bottom_text
                    horizontalalignment='center',
                    verticalalignment='top',
                    color=self.color_bottom_text)

        if show_IDs:
            for ind, level in enumerate(data):
                start = level[1]*(self.dimension+self.space)
                ax.text(start,
                        level[0]+self.offset,
                        str(ind),
                        horizontalalignment='right',
                        color='red')

        for idx, arrow in enumerate(self.arrows):
            for i in arrow:
                start = self.positions[idx]*(self.dimension+self.space)
                x1 = start + 0.5*self.dimension
                x2 = start + 0.5*self.dimension
                y1 = self.energies[idx]
                y2 = self.energies[i]
                gap = y1-y2
                gapnew = '{0:.2f}'.format(gap)
                middle= y1-0.5*gap
                ax.annotate("", xy=(x1,y1), xytext=(x2,middle), arrowprops=dict(color='green', width=2.5, headwidth=5))
                ax.annotate(s= gapnew, xy=(x2, y2), xytext=(x1, middle), color='green', arrowprops=dict(width=2.5, headwidth=5, color='green'), bbox=dict(boxstyle='round', fc='white'), ha='center', va = 'center')

        for idx, link in enumerate(self.links):
            for i in link:
                start = self.positions[idx]*(self.dimension+self.space)
                x1 = start + self.dimension
                x2 = self.positions[i[0]]*(self.dimension+self.space)
                y1 = self.energies[idx]
                y2 = self.energies[i[0]]

                gap = y2 - y1

                if i[5] == True:

                    x_middle = x1 + (x2-x1)/2
                    y_middle = y1 + gap / 2

                else:

                    x_middle = x1 + (x2-x1)/2
                    y_middle = y1 + gap / 2


                coef = '{0:.0f}'.format(i[4]*100)
                # linewidth = i[2]*abs(float(coef)/60)


                if i[4] > 0.05:
                    linewidth = i[2]*abs(float(coef)/60)
                else:
                    linewidth = i[2]*abs(float(coef)/20)

                ax.annotate("", xy=(x1,y1), xytext=(x_middle,y_middle), arrowprops={"arrowstyle":"-", "linestyle":i[1], "linewidth":linewidth, "color":i[3]})
                ax.annotate(text= coef+'%', xy=(x2, y2), xytext=(x_middle, y_middle), color=i[3], arrowprops={"arrowstyle":"-", "linestyle":i[1], "linewidth":linewidth,"color":i[3]}, ha='center', va = 'center')

        for box in self.electons_boxes:
            x, y, boxes, electrons, side, spacing_f = box
            plot_orbital_boxes(ax, x, y, boxes, electrons, side, spacing_f)

        self.ax = ax
        self.fig = fig

    def __auto_adjust(self,Energy_variation):
        # Max range between the energy
        # Energy_variation = abs(max(self.energies) - min(self.energies))
        if self.dimension == 'auto' or self.space == 'auto':
            # Unique positions of the levels
            # unique_positions = float(len(set(self.positions)))
            unique_positions = 3


            space_for_level = Energy_variation*self.ratio/unique_positions
            self.dimension = space_for_level*0.2
            self.space = space_for_level*0.8

        if self.offset == 'auto':
            self.offset = Energy_variation*self.offset_ratio
