class Category:
    def __init__(self, intervals):
        self.intervals=intervals

    def getIntervalCount(self, inPercentage=False, denominatorForEachInterval=None):
        """
        :param inPercentage: if we should show the values in percentage or not
        :param denominatorForEachInterval: To show in percentage, by default, we divide each interval count by the total count of all intervals. But you can also pass a vector with the denominator to use in each interval count
        :return:
        """
        if not inPercentage:
            return [ len(interval["data"]) for interval in self.intervals ]
        else:
            if denominatorForEachInterval==None:
                total = [sum( [ len(interval["data"]) for interval in self.intervals ] )] * len(self.intervals)
            else:
                total = denominatorForEachInterval
            return [float(len(interval["data"])) / float(total[index]) * 100.0 if total[index]>0 else 0 for index, interval in enumerate(self.intervals)]

    def __repr__(self):
        return self.intervals.__repr__()

    def __str__(self):
        return str(self.intervals)

    def __len__(self):
        return len(self.intervals)


class NumberCategory(Category):
    """
    Class that represents categories, as number intervals
    """
    def __init__(self, start, end, step):
        if type(start) is float or type(end) is float or type(step) is float:
            raise Exception("Do not user float as type to class Category, it is error-prone. Use int for integer intervals or decimal.Decimal for real intervals.")

        intervals=[]
        while start<end:
            intervals.append({"min": start, "max": start+step, "data": []})
            start+=step
        Category.__init__(self, intervals)

    def getCategoriesAsString(self, displayInterval=False, displayPlusOnFirstItem=False, displayPlusOnLastItem=False, displaySignal=True):
        """
        transform the intervals list into a list of string that will be the xlabels of the plot
        Basically, transform a vector like [-2, -1, 0, 1, 2] into ["-2+", "-1", "0", "+1", "+2+"]
        :return: list of string
        """
        intervalsAsString = []
        for i, interval in enumerate(self.intervals):
            if not displayInterval:
                lowerBound = interval["min"]
                # set the prefix
                prefix = "("
                if lowerBound > 0 and displaySignal:
                    prefix += "+"

                # set the suffix
                suffix = ""
                if (displayPlusOnFirstItem and i == 0) or (displayPlusOnLastItem and i == len(self.intervals) - 1):
                    suffix += "+"
                suffix += ")"

                intervalsAsString.append(prefix + str(lowerBound) + suffix)
            else:
                if i < len(self.intervals) - 1 or not displayPlusOnLastItem:
                    intervalsAsString.append("[%s,%s)" % (interval["min"], interval["max"]))
                else:
                    intervalsAsString.append("%s+"%(interval["min"]))

        return intervalsAsString

    def addDataPoint(self, point, dataToSave=None):
        """
        :param point: it is the value of the dataToSave - we use this to infer the category
        :return:
        """
        nbOfCategoriesItFit = 0

        if point < self.intervals[0]["min"]:
            self.intervals[0]["data"].append(dataToSave)
            nbOfCategoriesItFit += 1
        if point >= self.intervals[-1]["max"]:
            self.intervals[-1]["data"].append(dataToSave)
            nbOfCategoriesItFit += 1

        for i, interval in enumerate(self.intervals):
            if point >= interval['min'] and point < interval['max']:
                interval["data"].append(dataToSave)
                nbOfCategoriesItFit += 1

        if nbOfCategoriesItFit != 1:
            raise Exception("ERROR: %d fit %d categories..." % (point, nbOfCategoriesItFit))


class TextCategory(Category):
    """
    Class that represent categories, but as text
    """
    def __init__(self, categories):
        Category.__init__(self, [{"category": category, "data": []} for category in categories])

    def getCategoriesAsString(self, displayInterval=False, displayPlusOnFirstItem=False, displayPlusOnLastItem=False):
        return [interval["category"] for interval in self.intervals]

    def addDataPointAndIAlreadyKnowTheCategory(self, category, dataToSave=None):
        categoriesAsString = self.getCategoriesAsString()
        if category not in categoriesAsString:
            raise Exception("ERROR: non-existing category: %s..." % category)
        else:
            self.intervals[categoriesAsString.index(category)]["data"].append(dataToSave)
