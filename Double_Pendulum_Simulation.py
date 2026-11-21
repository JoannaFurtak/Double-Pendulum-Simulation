import math
import matplotlib.pyplot as plt
import numpy as np
import imageio.v2 as imageio
    
# parameters 
L1 = 7
L2 = 7
m1 = 7
m2 = 7
L = L1 + L2 + 0.1 #to na potrzebę wykresu, (granice wykresu)
T = 0
g = 5.5
theta01 = 0 #math.pi/1
alfa01 = 0
theta02 = math.pi/4
alfa02 = 0

t_max = 10
h = 0.01
N = int(t_max/h)
t = [0]

energyListRK4 = [0]
energyListEuler = [0]

# liczenie pochodnej z thety, wzory z m1 m2, l1, l2
def RHS(theta1, alfa1, theta2, alfa2):
    d1theta1 = alfa1
    d1theta2 = alfa2
    d2theta1 = float(-g*(2*m1+m2)*math.sin(theta1)- m2*g*math.sin(theta1-2*theta2)-2*m2*math.sin(theta1-theta2)*(d1theta2**2*L2+d1theta1**2*L1*math.cos(theta1 - theta2)))/(L1*(2*m1+m2-m2*math.cos(2*theta1-2*theta2)))
    d2theta2 = float(2*math.sin(theta1-theta2)*(d1theta1**2*L1*(m1+m2)+g*(m1+m2)*math.cos(theta1)+d1theta2**2*L2*m2*math.cos(theta1-theta2)))/(L1*(2*m1+m2-m2*math.cos(2*theta1 - 2*theta2)))
    return d1theta1, d2theta1, d1theta2, d2theta2

def energy_fun(alfa1, alfa2, theta1, theta2):
    y1 = -L1 * math.cos(theta1)
    y2 = y1 - L2 * math.cos(theta2)
    dx2 = alfa1*L1*math.cos(theta1)+alfa2*L2*math.sin(theta2)
    dy2 = alfa1*L1*math.sin(theta1)+alfa2*L2*math.sin(theta2)
    E = m1*(alfa1*L1)**2/2 + g*y1*m1 + m2*(dx2**2+dy2**2)/2 + g*y2*m2
    return E

X1 = [[theta01, alfa01, theta02, alfa02]] # lista rk4
X2 = [[theta01, alfa01, theta02, alfa02]] # lista euler

def euler(list2, i):
    ta = list2[i][0]
    ad = list2[i][1]
    th = list2[i][2]
    af = list2[i][3]
    k1 = h * np.array(RHS(ta, ad, th, af))
    listNext = list2[i] + k1
    list2.append(listNext.tolist())
    return list2

def RK4(list, i):
    ta = list[i][0]
    ad = list[i][1]
    th = list[i][2]
    af = list[i][3]

    k1 = h * np.array(RHS(ta, ad, th, af))
    k2 = h * np.array(RHS(ta + k1[0] / 2, ad + k1[1] / 2, th + k1[2] / 2, af + k1[3] / 2))
    k3 = h * np.array(RHS(ta + k2[0] / 2, ad + k2[1] / 2, th + k2[2] / 2, af + k2[3] / 2))
    k4 = h * np.array(RHS(ta + k3[0], ad + k3[1], th + k3[2], af + k3[3]))

    list_next = list[i] + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    list.append(list_next.tolist())
    return list

# rk4
x1 = [X1[0][0]] # d1theta1
y1 = [X1[0][1]] # d2theta1 
x2 = [X1[0][2]] # d1theta2
y2 = [X1[0][3]] # d2theta2

# euler
x1e = [X2[0][0]]
y1e = [X2[0][1]]
x2e = [X2[0][2]]
y2e = [X2[0][3]]

for i in range(N):
    X1 = RK4(X1, i)
    X2 = euler(X2, i)
    t.append(t[i]+h)
    x1.append(X1[i+1][0])
    y1.append(X1[i+1][1])
    x2.append(X1[i+1][2])
    y2.append(X1[i+1][3])
    x1e.append(X2[i+1][0])
    y1e.append(X2[i+1][1])
    x2e.append(X2[i+1][2])
    y2e.append(X2[i+1][3])

    energy = energy_fun(X1[i+1][1], X1[i+1][3], X1[i+1][0], X1[i+1][2])
    energyListRK4.append(energy)
    energy = energy_fun(X2[i+1][1], X2[i+1][3], X2[i+1][0], X2[i+1][2])
    energyListEuler.append(energy)
    
plt.plot(t, x1, label = " pierwsze ramię", color="purple") #pierwsze ramie
plt.plot(t, x2, label = "drugie ramię") 
plt.grid(True)
#plt.title("Wykres wartości przyspieszenia od czasu wyznaczony metodą RK4")
plt.ylabel("Odchylenie [pi*rad]")
plt.xlabel("Czas [s]")
plt.legend()
plt.show()

plt.plot(t, x1e, label = "pierwsze ramię", color="purple") # pierwsze ramie
plt.plot(t, x2e, label = "drugie ramię" )#color="green")
#plt.title("Wykres wartości przyspieszenia od czasu wyznaczony metodą Eulera")
plt.ylabel("Odchylenie [pi*rad]")
plt.xlabel("Czas [s]")
plt.grid(True) # FONT=
plt.legend(loc="lower left")
plt.show()

plt.plot(t, energyListRK4, label="RK4", color="red")
plt.plot(t, energyListEuler, label="Euler")
plt.ylabel("Energia [J]")
plt.xlabel("Czas [s]")
#plt.title("Wykres energii od czasu")
plt.grid(True)
plt.legend()
#plt.show()

# gif wahadlo podwojne
values_per_frame = 10 
t_ = t 
th1_ = x1 
th2_ = x2 

#przerobienie kątów na położenia w układzie kartezjańskim
x1_ = [L1 * math.sin(th1_[i]) for i in range(len(t_))]
y1_ = [-L1 * math.cos(th1_[i]) for i in range(len(t_))]
x2_ = [L1 * math.sin(th1_[i]) + L2 * math.sin(th2_[i]) for i in range(len(t_))]
y2_ = [-L1 * math.cos(th1_[i]) - L2 * math.cos(th2_[i]) for i in range(len(t_))]

frames = []
for i in range(len(t_)):
    if i % values_per_frame == 0:
        fig, ax = plt.subplots()
        ax.plot([0, x1_[i], x2_[i]], [0, y1_[i], y2_[i]], color='black') #rysowanie ramion wahadła
        ax.plot([x1_[i], x2_[i]], [y1_[i], y2_[i]], 'ko', color="blue") #rysowanie mas na wahadle
        ax.set_xlim(-L, L)
        ax.set_ylim(-L, L)
        ax.set_title(f'Wahadło o parametrach m1={m1}, m2={m2}, L1={L1}, L2={L2}')
        ax.text(-L, -L, f't = {t_[i]} sec') #podanie na wykresie upływającego czasu
        ax.axis("off") #wyłączenie tła wykresu (układu współrzędnych)
        filename = f'frame_{i+1}.png'
        plt.savefig(filename)
        frames.append(imageio.imread(filename))
        plt.close(fig)
    if i % 100 == 0 and i != 0: #śledzenie postępu programu
        print('progress =', i/len(t_))

imageio.mimsave('pendulum.gif', frames, duration=h*values_per_frame)
