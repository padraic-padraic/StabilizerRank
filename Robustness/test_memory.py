import resource
import threading
import time

from Utils.n_qubit_ops import gen_stabiliser_groups, gen_stabiliser_groups_tuple

class StoppableThread(threading.Thread):
    def __init__(self):
        super(StoppableThread, self).__init__()
        self.daemon = True
        self.__monitor = threading.Event()
        self.__monitor.set()
        self.__has_shutdown = False

    def run(self):
        '''Overloads the threading.Thread.run'''
        # Call the User's Startup functions
        self.startup()

        # Loop until the thread is stopped
        while self.isRunning():
            self.mainloop()

        # Clean up
        self.cleanup()

        # Flag to the outside world that the thread has exited
        # AND that the cleanup is complete
        self.__has_shutdown = True

    def stop(self):
        self.__monitor.clear()

    def isRunning(self):
        return self.__monitor.isSet()

    def isShutdown(self):
        return self.__has_shutdown

class MyLibrarySniffingClass(StoppableThread):
    def __init__(self, target_lib_call, arg):
        super(MyLibrarySniffingClass, self).__init__()
        self.target_function = target_lib_call
        self.arg = arg
        self.results = None

    def startup(self):
        # Overload the startup function
        print("Calling the Target Library Function...")

    def cleanup(self):
        # Overload the cleanup function
        print("Library Call Complete")

    def mainloop(self):
        # Start the library Call
        self.results = self.target_function(self.arg)

        # Kill the thread when complete
        self.stop()

if __name__ == "__main__":
    # Lib Testing Code
    mythread = MyLibrarySniffingClass(gen_stabiliser_groups_tuple, 3)
    mythread.start()

    start_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    delta_mem = 0
    max_memory = 0
    memory_usage_refresh = .005 # Seconds

    while(1):
        time.sleep(memory_usage_refresh)
        delta_mem = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) - start_mem
        if delta_mem > max_memory:
            max_memory = delta_mem

        # Uncomment this line to see the memory usuage during run-time 
        # print("Memory Usage During Call: %d MB" % (delta_mem / 1000000.0))

        # Check to see if the library call is complete
        if mythread.isShutdown():
            print('Done done done')
            # print(mythread.results)
            break;

    print("\nMAX Memory Usage in MB: " + str(round(max_memory / 1000.0, 3)))





