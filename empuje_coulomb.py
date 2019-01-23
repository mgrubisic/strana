from numpy import deg2rad, tan, sin, cos, rad2deg, pi, arctan, arange
import matplotlib.pyplot as plt
import seaborn as sns

def cohesionless_active(height, friction_angle_soil, unit_weight, wall_angle = 0, friction_angle_wall=0, angle_range=arange(40,85,1)):
    """no inclination, no seism, equivalent to mononobe
    enter angles in degrees"""
    alpha = deg2rad(wall_angle + friction_angle_wall)
    active_push = []
    for angle in angle_range:
        weight_soil = 0.5*height**2*unit_weight*(1.0/tan(deg2rad(angle)) + tan(deg2rad(wall_angle)))
        beta = deg2rad(angle - friction_angle_soil)
        active_push.append(weight_soil*sin(beta)/cos(alpha - beta))
    max_push = max(active_push)
    critical_angle = active_push.index(max_push) + angle_range[0]
    plt.figure()
    plt.minorticks_on()
    plt.title('Empuje activo sobre un muro de {}m de altura'.format(height))
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$E_a$ (kN)')
    # plt.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    # plt.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    # plt.grid(b=True, which='major', color='w', linewidth=2.0)
    # plt.grid(b=True, which='minor', color='w', linewidth=0.25)
    # plt.annotate(r'$\phi = {:.0f}^\circ$'.format(friction_angle_soil) + '\n' +  r'$\delta = {:.0f}^\circ$'.format(friction_angle_wall) +  '\n' +  r'$\lambda ={:.0f}^\circ $'.format(wall_angle), xytext = (0.8,0.5), textcoords = 'axes fraction')
    plt.text(critical_angle + 7, active_push[-1], r'$\phi = {:.0f}^\circ$'.format(friction_angle_soil) + '\n' +  r'$\delta = {:.0f}^\circ$'.format(friction_angle_wall) +  '\n' +  r'$\lambda ={:.0f}^\circ $'.format(wall_angle) + '\n\n' + r'$E_a = {:.1f}$ kN'.format(max_push) + '\n' + r'$\theta = {:.0f}^\circ$'.format(critical_angle))
    plt.plot(angle_range, active_push)
    plt.axhline(max_push, linestyle='--', color='k', linewidth = 0.1) # horizontal lines
    plt.axvline(critical_angle, linestyle='--', color='k', linewidth = 0.1) # vertical lines
    plt.savefig('talud-{}m-phi{}-delta{}-empuje-activo.pdf'.format(height, friction_angle_soil, friction_angle_wall))

    return max_push, critical_angle

def cohesionless_passive(height, friction_angle_soil, unit_weight, wall_angle = 0, friction_angle_wall=0, angle_range=arange(1,40)):
    """no inclination, no seism, equivalent to mononobe
    enter angles in degrees"""
    alpha = deg2rad(wall_angle - friction_angle_wall)
    print 'alpha = {}'.format(alpha)
    passive_push = []
    for angle in angle_range:
        weight_soil = 0.5*height**2*unit_weight*(1.0/tan(deg2rad(angle)) + tan(deg2rad(wall_angle)))
        beta = deg2rad(angle + friction_angle_soil)
        # passive_push.append(weight_soil*tan(beta)/(cos(alpha)+tan(beta)*sin(alpha)))
        passive_push.append( weight_soil*sin(beta)/cos(alpha - beta))
        print 'beta - {}'.format(beta)
        print sin(beta), cos(alpha-beta)
    max_push = min(passive_push)
    critical_angle = passive_push.index(max_push) + angle_range[0]
    plt.figure()
    plt.minorticks_on()
    plt.title('Empuje pasivo sobre un muro de {}m de altura'.format(height))
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$E_p$ (kN)')
    plt.text(critical_angle + 2, passive_push[-20], r'$\phi = {:.0f}^\circ$'.format(friction_angle_soil) + '\n' +  r'$\delta = {:.0f}^\circ$'.format(friction_angle_wall) +  '\n' +  r'$\lambda ={:.0f}^\circ $'.format(wall_angle) + '\n\n' + r'$E_p = {:.1f}$ kN'.format(max_push) + '\n' + r'$\theta = {:.0f}^\circ$'.format(critical_angle))
    plt.plot(angle_range, passive_push)
    plt.axhline(max_push, linestyle='--', color='k', linewidth = 0.1) # horizontal lines
    plt.axvline(critical_angle, linestyle='--', color='k', linewidth = 0.1) # vertical lines
    plt.savefig('talud-{}m-phi{}-delta{}-empuje-pasivo.pdf'.format(height, friction_angle_soil, friction_angle_wall))

    return max_push, critical_angle



def mononobe(height, friction_angle_soil, unit_weight, wall_angle = 0, friction_angle_wall=0, vertical_seismic_coeff = 0, horizontal_seismic_coeff = 0, fill_inclination = 0):
    phi = deg2rad(friction_angle_soil)
    alpha = deg2rad(fill_inclination)
    lambd = deg2rad(wall_angle)
    delta = deg2rad(friction_angle_wall)
    dzeta = arctan(horizontal_seismic_coeff/(1 - vertical_seismic_coeff))
    active_coeff = cos(phi - dzeta - lambd)**2/cos(dzeta)/cos(lambd)**2/cos(delta+lambd+dzeta)/(1+(sin(phi + delta)*sin(phi - dzeta-alpha)/cos(delta + lambd + dzeta)/cos(alpha - lambd))**.5)**2
    active_push = 0.5*height**2*unit_weight*(1 - vertical_seismic_coeff)*active_coeff

    return active_push


print cohesionless_active(6.0, 38.0, 18.9, 5, 0)
print cohesionless_active(6.0, 38.0, 18.9, 5, 10)
print cohesionless_active(6.0, 38.0, 18.9, 5, 20)

print cohesionless_passive(6.0, 38.0, 18.9, 5, 0)
print cohesionless_passive(6.0, 38.0, 18.9, 5, 10)
print cohesionless_passive(6.0, 38.0, 18.9, 5, 20)
