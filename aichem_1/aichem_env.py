import math
import seeding
import numpy as np
from qmwork_dftb import ChemEngine
# from plams_dftb3 import ChemEngine
# from plams_xtb import ChemEngine

bondlengthR1 = 2.992

class Env(object):

    # Set this in SOME subclasses
    metadata = {'render.modes': []}
    reward_range = (-float('inf'), float('inf'))
    spec = None

    # Set these in ALL subclasses
    action_space = None
    observation_space = None

    def step(self, action):

        raise NotImplementedError

    def reset(self):

        raise NotImplementedError

    def render(self, mode='human'):

        raise NotImplementedError

    def close(self):

        pass

    def seed(self, seed=None):

        return

    @property
    def unwrapped(self):
        """Completely unwrap this env.

        Returns:
            gym.Env: The base non-wrapped gym.Env instance
        """
        return self

    def __str__(self):
        if self.spec is None:
            return '<{} instance>'.format(type(self).__name__)
        else:
            return '<{}<{}>>'.format(type(self).__name__, self.spec.id)

    def __enter__(self):
        """Support with-statement for the environment. """
        return self

    def __exit__(self, *args):
        """Support with-statement for the environment. """
        self.close()
        # propagate exception
        return False

class ChemEnv(Env):
    """
    Description:
        A pole is attached by an un-actuated joint to a cart, which moves along a frictionless track. The pendulum starts upright, and the goal is to prevent it from falling over by increasing and reducing the cart's velocity.

    Source:
        This environment corresponds to the version of the cart-pole problem described by Barto, Sutton, and Anderson

    Observation:
        Type: Box(4)
        Num	Observation                 Min         Max
        0	Bondlength to R             -4.8            4.8
        1	Bondlength to P             -Inf            Inf
        2	Energy to R                -24 deg        24 deg
        3	Energy to p                 -Inf            Inf

    Actions:
        Type: Discrete(2)
        Num	Action
        0	Increase bondlength by 0.1 Amstrong
        1	Decrease bondlength by 0.1 Amstrong
        To: add another action to decrease the rate (0.1–0.01)

    Reward:
        Reward is 1 for every step taken, including the termination step, if energy is higher than previous step.

    Starting State:
        All observations are assigned a uniform random value in [bondlength of R –– bondlength of P]

    Episode Termination:
        Bondlength to R < 0
        Bondlength to P < 0
        Energy to R < 0
        Energy to P < 0
        Solved Requirements
        Considered solved when the average reward is greater than or equal to 195.0 over 100 consecutive trials.
    """

    def __init__(self):
        self.engine  = ChemEngine()
        self.bondChangerate = -0.1

        self.energy_threshold = 0
        self.bondlength_threshold = 0

        self.seed()
        self.state = None
        self.bondlength = None
        self.steps_beyond_done  = None
        self.previous_energy_r1 = None

    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]

    def step(self, action, step_name):

        bondChange = self.bondChangerate if action==1 else -self.bondChangerate
        self.bondlength +=  bondChange

        self.state = self.engine.run(self.bondlength, step_name)
        bondlength_r1, bondlength_r2, energy_r1, energy_r2 = self.state
        print("############################################################")
        print (step_name)
        print ('bondlength_r1, bondlength_r2, energy_r1, energy_r2')
        print (bondlength_r1, bondlength_r2, energy_r1, energy_r2)
        print ('self.previous_energy_r1')
        print (self.previous_energy_r1)

        done =  bondlength_r1 < self.bondlength_threshold \
                or bondlength_r2 < self.bondlength_threshold \
                or energy_r1 < self.energy_threshold \
                or energy_r2 < self.energy_threshold \
                or energy_r1 < self.previous_energy_r1
        done = bool(done)
        self.previous_energy_r1 = energy_r1
        if not done:
            reward = 1.0

        elif self.steps_beyond_done is None:
            # Pole just fell!
            steps_beyond_done = 0
            reward = 1.0
        else:
            if self.steps_beyond_done == 0:
                print ("You are calling 'step()' even though this environment has already returned done = True. You should always call 'reset()' once you receive 'done = True' -- any further steps are undefined behavior.")
            self.steps_beyond_done += 1
            reward = 0.0

        return np.array(self.state), reward, done, {}

    def reset(self):
        self.bondlength  =  bondlengthR1
        self.previous_energy_r1 = 0
        self.state = self.np_random.uniform(low=0.00, high=0.10, size=(4,))
        self.steps_beyond_done = None
        return np.array(self.state)