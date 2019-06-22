from orpytal import frames


class Trajectory(object):
    def __init__(self, states, name=None):
        self.states = states
        self.t_range = [self.states[0].t_since_rp, self.states[-1].t_since_rp]
        self.times = [st.t_since_rp for st in self.states]
        self.metadata = {}
        self.name = name
        self.frame = states[0].position.frame
        self.central_body = states[0].orbit.central_body

    def x_vals(self, frame='perifocal'):
        x_vals = []
        for st in self.states:
            frame_fn = getattr(st.position, frame)
            x_vals.append(frame_fn()[0])

        return x_vals

    def vals(self, frame_fn):
        x_vals = []
        y_vals = []
        z_vals = []

        if issubclass(frame_fn, frames.CoordinateFrame):
            frame_fn_name = frame_fn.fn_name
        else:
            frame_fn_name = frame_fn

        for st in self.states:
            frame_fn = getattr(st.position, frame_fn_name)
            pos = frame_fn(st).value
            x_vals.append(pos[0].m)
            y_vals.append(pos[1].m)
            z_vals.append(pos[2].m)

        return x_vals, y_vals, z_vals

    def inertial(self):
        return self.in_frame(frames.InertialFrame)

    def perifocal(self):
        return self.in_frame(frames.PerifocalFrame)

    def in_frame(self, frame_name):
        self.frame = frame_name

        self.x_vals, self.y_vals, self.z_vals = self.vals(frame_name)
        return self.x_vals, self.y_vals, self.z_vals

    def to(self, frame_name):
        return self.in_frame(frame_name)

    def start(self):
        self.states[0].set_vars()
        return self.states[0]

    def end(self):
        self.states[-1].set_vars()
        return self.states[-1]
