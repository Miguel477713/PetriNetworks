import os, numpy
from typing import Optional

if os.name == 'nt':
    os.system('cls')  # Clear the console for Windows
elif os.name == 'posix':
    os.system('clear')  # Clear the console for Linux/Mac

def MarkingKey(marking):
    """Return a hashable key representing the start marking."""
    # Convert numpy array to tuple of ints
    if marking is not None:
        return tuple(int(x) for x in marking)

class PetriNet:
    def __init__(self, pre, post, startMarking) -> None:

        if(len(pre) != len(post) or pre.shape[1] != post.shape[1]):
            raise ValueError("Pre and Post sizes are different")
        
        self.pre = pre
        self.post = post

        self.adjacencyMatrix = self.post - self.pre

        #print(self.adjacencyMatrix)

        self.startMarking = startMarking
        self.availableTransitions = self.GetAvailableTransistions()

    def MarkingKey(self):
        """Return a hashable key representing the start marking."""
        # Convert numpy array to tuple of ints
        return tuple(int(x) for x in self.startMarking)

    def GetAvailableTransistions(self) -> numpy.ndarray:
        availableTransistions = numpy.empty(self.pre.shape[1])

        for transition in range(self.pre.shape[1]):
            transitionIsAvailable = True
            for place in range(len(self.pre)):
                if self.startMarking[place] < self.pre[place][transition]:
                    transitionIsAvailable = False
                    break
            availableTransistions[transition] = int(transitionIsAvailable)
        return availableTransistions
    
    def GetNewMarking(self, firingTransition) -> Optional[numpy.ndarray]:
        if(firingTransition < 1 or firingTransition > self.pre.shape[1]):
            print(f"Transition t{firingTransition} is not defined")
            return None
        
        newMarking = numpy.empty(len(self.pre))
        firingTransition -= 1
        firingTransitionVector = numpy.zeros(self.pre.shape[1])
        firingTransitionVector[firingTransition] = 1

        if self.availableTransitions[firingTransition] == 1:
            newMarking = self.startMarking.T + self.adjacencyMatrix @ firingTransitionVector.T
        else:
            print(f"Transition t{firingTransition + 1} is not available")

        return newMarking
    
    def PrintAvailableTransitions(self) -> None:
        print(f"Available transitions from marking {self.startMarking} are: ")

        for i in range(len(self.availableTransitions)):
            if self.availableTransitions[i] == 1:
                print(f"t{i + 1}")


def main() -> None:
    pre = numpy.array([[1, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0],
                    [0, 0, 1, 0, 0],
                    [1, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 1],
                    ])

    post = numpy.array([[0, 0, 1, 0, 0],
                        [1, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0],
                        [0, 0, 0, 0, 1],
                        [1, 0, 0, 0, 0],
                        [0, 0, 0, 1, 0],
                    ])

    currentMarking = startMarking = numpy.array([1, 0, 0, 1, 0, 0])
    net = PetriNet(pre=pre, post=post, startMarking=startMarking)

    petriNetworks = {} #empty dict
    petriNetworks[net.MarkingKey()] = net


    userInput = 100
    while(True):
        net.PrintAvailableTransitions()
        userInput = input(f"Enter a transition [1..{pre.shape[1]}] or enter * to stop: ")

        if userInput == '*':
            break

        # try to parse transition number
        try:
            number = int(userInput)
        except ValueError:
            print("Invalid input, enter a number or * to stop.")
            continue

        currentMarking = net.GetNewMarking(number)
        if currentMarking is None:
            continue

        key = MarkingKey(currentMarking)
        if key in petriNetworks:
            net = petriNetworks[key]
        else:
            # otherwise create and store a new PetriNet for this marking
            net = PetriNet(pre=pre, post=post, startMarking=currentMarking)
            petriNetworks[net.MarkingKey()] = net

if __name__ == "__main__":
    main()
