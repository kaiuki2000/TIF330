import matplotlib.pyplot as plt
import numpy as np
plt.style.use('fast')

# After talking with the professor, I think this is the best way of doing this.
def SqrAbs(Psi):
    absPsi = []
    for i in Psi:
        absPsi.append(i[0]**2 + i[1]**2)
    return absPsi

# Definitions/Reading results
wVec       = np.genfromtxt(r'wVec.txt', delimiter=',')
# Psi        = np.genfromtxt(r'PsiXY.txt', delimiter=',')
PsiFT      = np.genfromtxt(r'PsiXY_FT.txt', delimiter=',')
AbsSqPsiFT = SqrAbs(PsiFT)

Spectrum = [AbsSqPsiFT[0]]
for i in range(1, int(len(wVec)/2.)):
    Spectrum.append(AbsSqPsiFT[i] + AbsSqPsiFT[len(wVec) - 1 - (i-1)])
# AbsSq.append(AbsSqPsiFT[int(len(wVec)/2.)]) # Ignoring the last \omega value, so our vector is symmetric.

plt.plot(wVec[:int(len(wVec)/2)], Spectrum, label = r"Spectrum: $|C_{\omega}|^2 + |C_{-\omega}|^2$")

# plt.plot(wVec[:int(len(wVec)/2)], [PsiFT[i][0] for i in range(0, len(PsiFT))][:int(len(wVec)/2)], label = "Real part")
# plt.plot(wVec[:int(len(wVec)/2)], [PsiFT[i][1] for i in range(0, len(PsiFT))][:int(len(wVec)/2)], label = "Imaginary part")
# plt.plot(wVec[int(len(wVec)/2):][::-1], [PsiFT[i][0] for i in range(0, len(PsiFT))][int(len(wVec)/2):][::-1], color = 'C0')
# plt.plot(wVec[int(len(wVec)/2):][::-1], [PsiFT[i][1] for i in range(0, len(PsiFT))][int(len(wVec)/2):][::-1], color = 'C1')

# Peak around 3; Colour map of Psi (Real part is enough).
# Absolute value plot
# plt.plot(wVec[:int(len(wVec)/2)], AbsSqPsiFT[:int(len(wVec)/2)], label = 'Absolute value')
# plt.plot(wVec[int(len(wVec)/2):], AbsSqPsiFT[int(len(wVec)/2):], color = 'C0')
# plt.xlim(-2.5, 5.)

y = [0., 2.e-5, 4.e-5, 6.e-5, 8.e-5, 10.e-5, 12.e-5, 14.e-5]
ylabels = [f'{label:,}' for label in y]
# plt.yticks(y, ylabels)

plt.xlabel(r"Frequency: $\omega$")
# plt.ylabel(r"$\mathscr{F}(|\psi(0.1, 0.0)|^2)$")
# plt.ylabel(r"$\mathscr{F}(\psi(0.1, 0.0))$")
plt.ylabel(r"$|C_{\omega}|^2 + |C_{-\omega}|^2$")
plt.title(r"Spectrum of $\psi(t)$ at $(x, y) = (0.1, 0.0)$")
plt.grid(True)
plt.legend()
# plt.xlim(0., 5.)
# plt.ylim(0., 1.)
plt.tight_layout()
plt.savefig('Spectrum.png', dpi = 500)

print(f"Fundamental frequency = wVec[{np.argmax(Spectrum)}] = {wVec[np.argmax(Spectrum)]}")

# Ask Teacher/Nico/Viktor exactly what I should plot here. Are my peaks good?
# They seem a bit 'weak'... (In the case of F{|Psi|^2}).
# Looking back, I think it makes more sense to plot the FT of the full Psi; Since, H\Psi = E\Psi.
# Eigenvalues are of \Psi, not |\Psi|^2.
# I should also try Problem 2 again.