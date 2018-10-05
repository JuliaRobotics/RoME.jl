# Direct unimodal bearing range only calculation

# ground truth
wX1 = [1;0;0] # [x,y,th]
wL1 = [-1;3]
wL2 = [3;3]
wL3 = [3;0]

# cartesian measurements
DX1 = wL1-wX1[1:2]
DX2 = wL2-wX1[1:2]
DX3 = wL3-wX1[1:2]

# find location
# wTb = [R(th) tr; 0 0 1]
# Use linear relationship: bDX = bTw*wL
bDX = ones(3,3)
bDX[1:2,1] = DX1
bDX[1:2,2] = DX2
bDX[1:2,3] = DX3

wL = ones(3,3)
wL[1:2,1] = wL1
wL[1:2,2] = wL2
wL[1:2,3] = wL3

bTw = (wL'\bDX')' # use transpose for pivoted Cholesky inverse
wTb = inv(bTw)

# answers
wThetab = atan(wTb[1,2], wTb[1,1])
wTransb = wTb[1:2,3]

# Brute force least squares: (A'A)\(A')*y
# https://en.wikipedia.org/wiki/Linear_least_squares_(mathematics)

wL4 = [0;-3]
DX4 = wL4-wX1[1:2]

bDX = ones(3,4)
bDX[1:2,1] = DX1
bDX[1:2,2] = DX2
bDX[1:2,3] = DX3
bDX[1:2,4] = DX4

wL = ones(3,4)
wL[1:2,1] = wL1
wL[1:2,2] = wL2
wL[1:2,3] = wL3
wL[1:2,4] = wL4

# bDX*bDX' = bDX*wL'*bTw'
bTw = ((bDX*wL')\(bDX*bDX'))'
