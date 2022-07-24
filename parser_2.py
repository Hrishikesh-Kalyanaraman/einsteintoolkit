from collections import defaultdict
import re
def create_summary(file):
    '''
        This function parses the test results from the build__2_1.log file
        into a dictionary
    '''
    sum_data = {}
    with open(file, "r") as fp:
        lines = fp.read().splitlines()
        line_number = 0
        lines_2 = lines[0].split(" ")
        sum_data["Total available tests"] = int(lines_2[1])
        sum_data["Unrunnable tests"] = 0
        sum_data["Runnable tests"] = 0
        sum_data["Number of tests passed"] = 0
        sum_data["Number failed"] = 0
        while line_number < len(lines) - 1:
            #print(line_number)
            while not "Running Test" in lines[line_number]:
                line_number += 1
                if (line_number == len(lines) - 1):
                    break
            if (line_number < len(lines) - 1):
                test_result = lines[line_number + 1].split(":")[0]
                if (test_result == "Passed"):
                    sum_data["Number of tests passed"] += 1
                    sum_data["Runnable tests"] += 1
                elif (test_result == "Unrunnable"):
                    sum_data["Unrunnable tests"] += 1
                else:
                    sum_data["Number failed"] += 1
                    sum_data["Runnable tests"] += 1
                line_number += 1
        #print(sum_data)
    return sum_data

#
def get_tests(readfile):
    '''
        This functions reads the names of the tests that failed and passed
        from the log files and returns the set of them in tuples.
    '''
    passed = set()
    failed = set()
    with open(readfile, "r") as fp:
        lines = fp.read().splitlines()
        line_number = 0
        while line_number < len(lines) - 1:
            #print(line_number)
            while not "Running Test" in lines[line_number]:
                line_number += 1
                if (line_number == len(lines) - 1):
                    break
            if (line_number < len(lines) - 1):
                test_result = lines[line_number + 1].split(":")[0]
                if (test_result == "Passed"):
                    passed.add(lines[line_number].split("/")[1])
                else:
                    failed.add(test_result)
                line_number += 1
        #print(passed, failed)
    return passed, failed

def test_comp(readfile_new, readfile_old):
    '''
        This function compares the tests from the previous run
        to see which failed, which failed this time but did
        not in the last, which passed this time but not in the last,
        and the new and removed tests.
    '''
    passed_n, failed_n = get_tests(readfile_new)
    passed_o, failed_o = get_tests(readfile_old)
    newly_p = passed_n - passed_o
    newly_f = failed_n - failed_o
    new_tests = (passed_n.union(failed_n)) - (passed_o.union(failed_o))
    missing_tests = (passed_o.union(failed_o)) - (passed_n.union(failed_n))
    types = ["Failed Tests", "Newly Passing Tests", "Newly Failing Tests", "Newly Added Tests", "Removed Tests"]
    tests = [failed_n, newly_p, newly_f, new_tests, missing_tests]
    test_dict = {types[i]: tests[i] for i in range(len(types))}
    return test_dict


def get_data(name):
    '''
        Retrieves singular field of data from the data csv as a list
    '''
    data = {}
    with open('new_test_nums.csv',encoding='utf-8-sig') as csvfile:
        fields = csvfile.readline().strip().split(",")
        #g(fields)
        name_i = fields.index(name)
        line = csvfile.readline()
        while line != "":
            entry = line.strip().split(",")
            data[entry[0]] = float(entry[name_i])
            line = csvfile.readline()
    return data