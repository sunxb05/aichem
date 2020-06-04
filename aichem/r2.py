import argparse
import numpy as np
import random
from aichem_env import ChemEnv
from ddqn_game_model import DDQNTrainer, DDQNSolver


class rl:

    def __init__(self):
        game_name, game_mode, total_step_limit, total_run_limit = self._args()
        env = ChemEnv()
        self.observation_space = 4
        self.action_space      = 20
        self.highenergy        = [3,3,3,3,3]
        # observation_space = env.observation_space.shape[0]
        # action_space = env.action_space.n
        self._main_loop(self._game_model(game_mode, game_name, self.observation_space, self.action_space), env, total_step_limit, total_run_limit)

    def _main_loop(self, game_model, env, total_step_limit, total_run_limit):
        run = 0
        total_step = 0
        while True:
            if total_run_limit is not None and run >= total_run_limit:
                print ("Reached total run limit of: " + str(total_run_limit))
                exit(0)

            run += 1
            state = env.reset()
            state = np.reshape(state, [1, self.observation_space])
            step = 0
            score = 0
            while True:
                if total_step >= total_step_limit:
                    print ("Reached total step limit of: " + str(total_step_limit))
                    exit(0)
                total_step += 1
                step += 1
                action = game_model.act(state)
                state_next_list, reward, terminal, info = env.step(action,step,sum(self.highenergy)/5)
                score += reward
                state_next = np.reshape(state_next_list, [1, self.observation_space])
                game_model.remember(state, action, reward, state_next, terminal)
                state = state_next
                game_model.step_update(total_step)
                if terminal:
                    self.highenergy.sort()
                    self.highenergy[0] = state_next_list[3]
                    game_model.save_run(score, step, run)
                    break

    def _args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-g", "--game", help="Choose from available games: aichem. Default is 'aichem'.", default="aichem")
        parser.add_argument("-m", "--mode", help="Choose from available modes: ddqn_train, ddqn_test, . Default is 'ddqn_training'.", default="ddqn_training")
        parser.add_argument("-tsl", "--total_step_limit", help="Choose how many total steps (frames visible by agent) should be performed. Default is '100000'.", default=1000000, type=int)
        parser.add_argument("-trl", "--total_run_limit", help="Choose after how many runs we should stop. Default is 10000 (no limit).", default=100000, type=int)
        args = parser.parse_args()
        game_name = args.game
        game_mode = args.mode
        total_step_limit = args.total_step_limit
        total_run_limit = args.total_run_limit

        print ("Selected game: " + str(game_name))
        print ("Selected mode: " + str(game_mode))
        print ("Total step limit: " + str(total_step_limit))
        print ("Total run limit: " + str(total_run_limit))
        return game_name, game_mode, total_step_limit, total_run_limit

    def _game_model(self, game_mode, game_name,observation_space, action_space):
        if game_mode == "ddqn_training":
            return DDQNTrainer(game_name, observation_space, action_space)
        elif game_mode == "ddqn_testing":
            return DDQNSolver(game_name, observation_space, action_space)
        else:
            print ("Unrecognized mode. Use --help")
            exit(1)


if __name__ == "__main__":
    rl()
