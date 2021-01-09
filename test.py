# -*-coding:utf-8-*-
# 2D_FDTD

import numpy
import math
import matplotlib.pyplot as plt

#初期設定値（定数）
dimension = 2            #次元
rhow = 1000.0          #水密度
frequency = 1.0e+6  #周波数 1Mhz なんでも
vertical_wave_sound_velocity=1500.0 #水の縦波音速
T = 1/frequency
rhob = 2000.0


#分解能　大きければ大きいほど精度が上がる
#設定値によってはPCがフリーズするので注意
dx = 1e-5
dy = dx
dt = dx/vertical_wave_sound_velocity/math.sqrt(dimension)

#シミュレーション空間の範囲
model = numpy.zeros((500, 500))
nn = model.shape
nx = nn[0]-1
ny = nn[1]-1
nt = 1000


# 入力波形 今回は、sin 1波
n = int(T / dt)
wave2 = numpy.array([math.sin(math.pi * 2 / T * dt * ms)
                     for ms in range(0, n+2)])

# 音圧と音速行列を定義
Pxx = numpy.zeros((nx, ny), dtype=float)
Ux = numpy.zeros((nx+1, ny), dtype=float)
Uy = numpy.zeros((nx, ny+1), dtype=float)

# 音波発生位置 今回は、真ん中
inputpointx = int(nx/2)
inputpointy = int(ny/2)

for s in range(nt):

    # 入力する波
    if s <= n+2:
        Pxx[inputpointx, inputpointy] = wave2[s-1]
    else:
        Pxx[inputpointx, inputpointy] = 0

    # 波の伝わり方の計算
    Ux[1: nx, 0: ny] = Ux[1: nx, 0: ny] - \
        (1/rhow) * (dt/dx) * (Pxx[1: nx, 0: ny] - Pxx[0: nx-1, 0: ny])
    Uy[0: nx, 1: ny] = Uy[0: nx, 1: ny] - \
        (1/rhow) * (dt/dy) * (Pxx[0: nx, 1: ny] - Pxx[0: nx, 0: ny-1])

    Pxx[0: nx, 0: ny] = Pxx[0: nx, 0: ny] - rhow * \
        vertical_wave_sound_velocity ** 2 * dt * ((1/dx)*(Ux[1: nx+1, 0: ny]-\
        Ux[0: nx, 0: ny]) + (1/dy)*(Uy[0: nx, 1: ny+1]-Uy[0: nx, 0: ny]))

#シミュレーションの表示
    if s % 5 == 0:
        print("count: {0}".format(s))
        plt.cla()
        plt.imshow(Pxx, vmin=-
               1, vmax=1, origin=0)
        plt.draw()
        plt.pause(0.00001)

