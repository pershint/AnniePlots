import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set_context("poster")

#First, range our theta from -45 to 45

t_low = 21.
t_high = 3.
t_c = 42.

def getConeCharge(theta_arr):
    theta_diff = np.array(theta_arr - t_c)
    low_ind = np.where(theta_diff > 0)[0]
    high_ind = np.where(theta_diff < 0)[0]
    #not fast, but works
    cone_charge = []
    for j,theta in enumerate(theta_diff):
        if j in high_ind:
            p = (0.75 + (0.25/(1.+(theta**2/t_low**2))))
            cone_charge.append(p)
        else:
            p = (1./(1.+(theta**2/t_high**2)))
            cone_charge.append(p)
    return cone_charge

if __name__=='__main__':
    theta = np.arange(0.,50.,0.5)
    cone_charge = getConeCharge(theta)
    sns.set_style("whitegrid")
    sns.axes_style("darkgrid")
    xkcd_colors = ['eggplant', 'grass', 'vomit', 'warm pink', 'green', 'grass']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(theta,cone_charge, linewidth=5)
    plt.title("Cone charge scaling factor",fontsize=34)
    for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(30)
    for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(30)
    plt.xlabel(r"$\delta_{i}$ (degrees)",fontsize=34)
    plt.show()
