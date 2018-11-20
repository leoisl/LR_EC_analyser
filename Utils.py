#Some utilities classes

class FileUtils:
    '''
    Some files for easy managing of files (pretty few stuff really)
    '''
    @staticmethod
    def readFileComposedOfPairStringIntToDict(filename):
        stats = {}
        with open(filename) as file:
            for line in file:
                lineSplit = line.split()
                stats[lineSplit[0]] = int(lineSplit[1])
        return stats

    @staticmethod
    def readFileComposedOfPairIntsToDict(filename):
        stats = {}
        with open(filename) as file:
            for line in file:
                lineSplit = line.split()
                stats[int(lineSplit[0])] = int(lineSplit[1])
        return stats