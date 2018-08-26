from orpytal import frames

class Trajectory(object):
    name = None

    def __init__(self, states, name=None):
        self.states = states
        self.t_range = [self.states[0].t_since_rp, self.states[-1].t_since_rp]
        self.times = [st.t_since_rp for st in self.states]
        self.metadata = {}
        self.name = name
        self.frame = states[0].position.frame

    def x_vals(self, frame='orbit_fixed'):
        x_vals = []
        for st in self.states:
            frame_fn = getattr(st.position, frame)
            x_vals.append(frame_fn()[0])

        return x_vals

    def vals(self, frame):
        x_vals = []
        y_vals = []
        z_vals = []

        for st in self.states:
            frame_fn = getattr(st.position, frame)
            pos = frame_fn().value
            x_vals.append(pos[0].m)
            y_vals.append(pos[1].m)
            z_vals.append(pos[2].m)

        return x_vals, y_vals, z_vals

    def inertial(self):
        self.frame = frames.InertialFrame
        return self.in_frame('inertial')

    def orbit_fixed(self):
        self.frame = frames.OrbitFixedFrame
        return self.in_frame('orbit_fixed')

    def in_frame(self, frame_name):
        self.x_vals, self.y_vals, self.z_vals = self.vals(frame_name)
        return self.x_vals, self.y_vals, self.z_vals

    def start(self):
        return self.states[0]

    def end(self):
        return self.states[-1]
