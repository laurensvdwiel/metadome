from joblib.parallel import cpu_count

def CalculateNumberOfActiveThreads(numberOfTasks):
    if(cpu_count() == 2):
        return cpu_count()
    elif numberOfTasks < cpu_count():
        return numberOfTasks
    else:
        return cpu_count()