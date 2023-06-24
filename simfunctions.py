import numpy as np
R = 1e-6
kb = 1.38*10**-23
viscosity = 1*10**-3
#D_T = 1 * 1e-13
D_T = kb*300/(6*np.pi*viscosity*R)
D_R = 1
v = 2e-6
v0 = 50 * 1e-6
delta_t = 0.03
variance = 1
parameters = np.array([R, D_T, D_R, v, delta_t, variance])
box_size = 100 * 1e-6
cut_off_distance = 8 * R

phoreticConstant = v0 * R ** 2
twoPi = 2*np.pi
twoR = 2*R
#Pre-compute certain constants for simulate_motion:
omega = np.sqrt(2*D_T*delta_t)
beta = np.sqrt(2*D_R*delta_t)


def initialize_positions(_position_list, _n_particles, _box_size):
    for i in range(_n_particles):
        _position_list[i, :] = np.random.rand(2) * _box_size
    return _position_list


def initialize_angles(_phi_list, _n_particles):
    for i in range(_n_particles):
        _phi_list[i] = np.random.rand() * 2 * np.pi
    return _phi_list


def simulate_motion(_position, _phi, _phoretic_velocity, _delta_t, omega ,beta, _v, _v0, _variance, _n_timesteps):
    _new_position = np.zeros(2)
    _W_phi = np.random.normal(loc=0, scale=variance)
    _W_x = np.random.normal(loc=0, scale=variance)
    _W_y = np.random.normal(loc=0, scale=variance)
    _angle = _phi + _W_phi * beta

    #Update X and Y respectively:
    _new_position[0] = _position[0] + _v * np.cos(_angle) * _delta_t + omega * \
                       _W_x + _phoretic_velocity[0] * delta_t

    _new_position[1] = _position[1] + _v * np.sin(_angle) * _delta_t + omega * \
                       _W_y  + _phoretic_velocity[1] * delta_t


    return np.array(_new_position)


def phoretic_interaction(_v0, _R, _position1, _position2, _cut_off_distance):
    _distance = np.linalg.norm(_position2 - _position1)
    
    if (_distance < _cut_off_distance) and (_distance != 0):
        _direction = (_position2 - _position1) / _distance
        _phoretic_velocity = (phoreticConstant /( _distance ** 2)) * _direction
        return _phoretic_velocity
    else:
        _phoretic_velocity = np.array([0.0, 0.0])
        return _phoretic_velocity


def hard_sphere_correction(_position1, _position2):
    _center_distance = np.linalg.norm(_position1 - _position2)
    _overlap = twoR - _center_distance
    if (_overlap > 0) and (_overlap != twoR):
        _distance_to_move = _overlap / 2
        _direction_to_move = (_position1 - _position2) / _center_distance

        
        _new_position1 = _position1[:] + _direction_to_move * _distance_to_move
        _new_position2 = _position2[:] - _direction_to_move * _distance_to_move
        return _new_position1, _new_position2
  
    else:
        return _position1, _position2



def check_overlap(_position_list, _R, _n_particles):
    for i in range(_n_particles):
        for k in range(i + 1, _n_particles):
            _position1, _position2 = hard_sphere_correction(_position_list[i, :], _position_list[k, :])
            _position_list[i, :] = _position1
            _position_list[k, :] = _position2
    return _position_list

# def check_overlap(_position_list,_R, _n_particles):
#     for i in range(_n_particles):
#         for k in range(_n_particles):
#             if i != k:
#                 _position1, _position2 = hard_sphere_correction(_position_list[i, :], _position_list[k, :])
#                 _position_list[i, :] = _position1
#                 _position_list[k, :] = _position2
#     return _position_list



def outside_box(_position_list, _n_particles, _box_size):
    _discontinuity_pos = np.ones(_n_particles)
    for i in range(_n_particles):
        if _position_list[i, 0] > _box_size:
            _position_list[i, 0] = _position_list[i, 0] - _box_size
            _discontinuity_pos[i] = np.nan
        elif _position_list[i, 0] < 0:
            _position_list[i, 0] = _position_list[i, 0] + _box_size
            _discontinuity_pos[i] = np.nan
        if _position_list[i, 1] > _box_size:
            _position_list[i, 1] = _position_list[i, 1] - _box_size
            _discontinuity_pos[i] = np.nan
        elif _position_list[i, 1] < 0:
            _position_list[i, 1] = _position_list[i, 1] + _box_size
            _discontinuity_pos[i] = np.nan
    return _position_list, _discontinuity_pos
