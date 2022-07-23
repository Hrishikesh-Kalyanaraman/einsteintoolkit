#######################################################################################################################
from requests.api import get
import sys
import os, csv, time, requests, math
from datetime import datetime
import matplotlib.pyplot as plt

# Different Bokeh modules
from bokeh.models.annotations import Legend
import bokeh.plotting as bplt
import bokeh.models.tools as btools
from bokeh.models import Panel, Tabs, DataRange1d
import bokeh.models.callbacks as bcall
from bokeh.resources import CDN
from bokeh.embed import components
from bokeh.layouts import row

from store_2 import get_version, copy_index, get_commit_id
from bokeh.palettes import viridis
from bokeh.transform import factor_cmap
from bokeh.resources import CDN
from bokeh.embed import file_html

# to generate a commit log
from pygit2 import Repository, Oid
from pygit2 import GIT_SORT_TOPOLOGICAL, GIT_SORT_REVERSE
from datetime import datetime, timezone, timedelta
import time

# Functions from parser.py
from parser_2 import create_summary, get_tests, test_comp, get_data
import glob

records = os.listdir("./records")
curr_ver = get_version() - 1
curr = f"./records/version_{curr_ver}/build__2_1_{curr_ver}.log"
last = f"./records/version_{curr_ver - 1}/build__2_1_{curr_ver - 1}.log"

repo = Repository('.git')
baseurl = repo.remotes["origin"].url.replace("git@", "https://").replace(".git", "")
#print(baseurl)


def gen_commits():
    '''
        This function generates a list of commits that have been made since the last run
        If the workflow was run manually it will say so.
    '''

    # TODO: turn into convenience function
    curr_commit_id = Oid(hex=get_commit_id(curr_ver))
    last_commit_id = Oid(hex=get_commit_id(curr_ver - 1))
    commits = []
    for commit in repo.walk(curr_commit_id, GIT_SORT_TOPOLOGICAL):
        if (commit.id == last_commit_id):
            break
        commits.append(commit)

    lines = []
    for count, commit in enumerate(commits):
        message = commit.message
        message = message.replace("\n\n", "\n")
        message = message.replace('\n', '<br>')
        tzinfo = timezone(timedelta(minutes=commit.author.offset))
        # print(tzinfo)
        dt = datetime.fromtimestamp(commit.author.time, tz=timezone.utc)
        # print(dt)
        date = str(dt)
        lines.append("<tr> <td> <b> Commit: </b> </td> <td> <b> " + str(count + 1) + " </b> </td> </tr> \n")
        lines.append("<tr> <td> Date: (In UTC)  </td> <td> " + date + "  </td> </tr> \n")
        lines.append("<tr> <td> Message: </td> <td>" + message + "</td> </tr> \n")
        # print(lines)
    return "\n".join(lines)


# log_link=f"{baseurl}/blob/master/records/version_{last_ver}/{last[ext-1:]+str(last_ver+1)}"
def gen_diffs(readfile):
    '''
        This function generates the html table that shows
        the comparison of test logs from last version generated
        by test_comp
    '''

    # The test_comp function provides tests that failed, were newly added or newly removed
    test_comparison = test_comp(readfile, last)
    if len(test_comparison["Failed Tests"]) != 0:
        print("TESTS_FAILED=True")

    # Setup the header for the table
    output = '''<table class="table table-bordered " >
    <caption style="text-align:center;font-weight: bold;caption-side:top">Failed Tests and Changes</caption>
    <tr><th></th><th>logs(1_process)</th><th>logs(2_processes)</th><th>diffs(1_process)</th><th>diffs(2_processes)</th></tr>\n'''

    for result in test_comparison.keys():
        # For each test make a header with the description of why that test is being shown(failed, newly added, newly failing)
        output += f'''<tr><th colspan="5">''' + result + "</th></tr>\n"

        # If no such test exists add empty row
        if (len(test_comparison[result]) == 0):
            output += "<tr><td></td></tr>"

        # For each test get the thorn name and the current version
        for test in test_comparison[result]:
            thorn = test.split()[-1]
            thorn = thorn[:len(thorn) - 1]
            test_name = test.split()[0]
            ver = curr_ver

            # Since the removed test would have been stored in the curr version subtract 1 from the version number
            if ("Removed" in result):
                ver -= 1

            # Links for logs and diffs of the tests in the test_comparison dictionary based on the number of procs
            # logl1 = f"{baseurl}/tree/gh-pages/records/version_{ver}/sim_{ver}_1/{thorn}/{test_name}.log"
            # logl2 = f"{baseurl}/tree/gh-pages/records/version_{ver}/sim_{ver}_2/{thorn}/{test_name}.log"
            # diffl1 = f"{baseurl}/tree/gh-pages/records/version_{ver}/sim_{ver}_1/{thorn}/{test_name}.diffs"
            # diffl2 = f"{baseurl}/tree/gh-pages/records/version_{ver}/sim_{ver}_2/{thorn}/{test_name}.diffs"

            # Check if these files are available if not display not avaible on the table
            # if (os.path.isfile("./" + logl1[logl1.find("records"):])):
            #     output += f"  <tr><td>{test}</td><td><a href='{logl1}'>log</a></td>"
            # else:
            #     output += f" <tr><td>{test}</td><td>Not Available</td>"
            # if (os.path.isfile("./" + logl2[logl2.find("records"):])):
            #     output += f"  <td><a href='{logl2}'>log</a></td>"
            # else:
            #     output += f" <td>Not Available</td>"
            # if (os.path.isfile("./" + diffl1[diffl1.find("records"):])):
            #     output += f"<td><a href='{diffl1}'>diff</a></td>"
            # else:
            #     output += f"<td>Not Available</td>"
            # if (os.path.isfile("./" + diffl2[diffl2.find("records"):])):
            #     output += f"<td><a href='{diffl2}'>diff</a></td></tr>\n"
            # else:
            #     output += f"<td>Not Available</td></tr>\n"

    output += "</table>"
    return output



def plot_test_data(readfile):
    # Get dataa from the csv and create lists for each field
    runnable = list(get_data("Runnable tests").values())
    times = list(get_data("Number of tests passed").keys())
    passed = list(get_data("Number of tests passed").values())
    #time_taken = list(get_data("Time Taken").values())
    #compile_warn = list((get_data("Compile Time Warnings").values()))
    build_no = list(get_data("Build Number").values())

    # Get the of dictionary of thorns with their warning counts
    #warning_thorns = get_warning_thorns(f"records/version_{curr_ver}/build_{curr_ver}.log")

    # Turn that dictionary into lists so you can pick the thorns with most warnings
    # counts = list(warning_thorns.values())
    # counts_trunc = sorted(counts, reverse=True)[:7]
    # warning_types_trunc = []
    # warning_types_list = list(warning_thorns.keys())
    #
    # for count in counts_trunc:
    #     i = counts.index(count)
    #     warning_types_trunc.append(warning_types_list[i])
    #     warning_types_list.pop(i)
    #     counts.pop(i)
    # counts = counts_trunc
    # warning_types_list = warning_types_trunc

    axis = []  # values displayed on the x axis

    Date = []  # dates in YYYY/MM/DD form
    for i in range(len(times)):
        date = datetime.fromtimestamp(float(float(times[i])), tz=None)  # changing to UNIX timestamp from float
        axis.append(date)
        Date.append(date.strftime("%Y-%m-%d"))
        # print(axis)
    print(axis)
    # The python library bokeh has a special data structure called a column data source that functions similarly
    # to a dictionary
    src = bplt.ColumnDataSource(data=dict(
        t=times,
        nt=axis,
        d=Date,
        b=build_no,
        rt=runnable,
        tp=passed,
        # timet=time_taken,
        # cmt=compile_warn,
        xax=[0] * len(times),
        #url=[f"./index_{x + 1}.html" for x in range(0, curr_ver)],
    ))

    # p is the first figure an area chart with the number of tests passed out of the ones ran
    # Tools attribute gives ways to manipulate the plot such as having clickable points, scrool to zoom and pan to zoom.
    # The rest of the attributes should be self explanatory
    # start = times[len(axis) - 10]
    # end = times[len(axis) - 1]
    # print(start, end)
    p = bplt.figure(x_range= (axis[len(axis) - 10], axis[len(axis) - 1]), y_range=(0, 12), plot_width=1000, plot_height=600
                    , tools="hover,tap,xwheel_zoom,box_zoom,reset",
                    y_axis_label="Number of Tests", x_axis_label="Date", x_axis_type="datetime",
                    title="Passed Tests", toolbar_location="below", sizing_mode='scale_width', tooltips=[("Date", "@d"),
                                                                                                         (
                                                                                                         "Runnable Tests",
                                                                                                         "@rt"), (
                                                                                                         "Passing Tests",
                                                                                                         "@tp"),
                                                                                                         ("Build Number",
                                                                                                          "@b")])
    # tooltips are used for the hover tool

    # Circles are points on the graph
    p.circle('nt', 'rt', size=10, color="green", source=src, legend_label="Runnable Tests")
    p.circle('nt', 'tp', size=10, color="blue", source=src, legend_label="Number of Tests Passed")

    # The taptool helps have these points link to the previous builds
    # url = "@url"
    # taptool = p.select(type=btools.TapTool)
    # taptool.callback = bcall.OpenURL(url=url)

    # This part fills in the area below the points
    p.varea(y1='rt', y2='xax', x='nt', color="green", source=src, alpha=0.5)
    p.varea(y1='tp', y2='xax', x='nt', color="blue", source=src, alpha=0.5)

    # The graphs are displayed in a tabs and this part sets that up
    tab1 = Panel(child=p, title="Test Results")
    p.legend.location = "top_left"

    # # This graph is for how long the testing part takes uses similar code to the first one but instead of area it has lines connecting points
    # p1 = bplt.figure(y_range=(0, max(time_taken) + 5), plot_width=1000, plot_height=600,
    #                  tools="tap,wheel_zoom,box_zoom,reset",
    #                  y_axis_label="Time(minutes)", x_axis_label="Date", x_axis_type="datetime",
    #                  title="Time Taken for Tests", toolbar_location="below", sizing_mode='scale_width')
    # p1.circle('nt', 'timet', size=10, color="blue", source=src)
    # p1.line('nt', 'timet', color="blue", source=src)
    # taptool = p1.select(type=btools.TapTool)
    # taptool.callback = bcall.OpenURL(url=url)
    # tab2 = Panel(child=p1, title="Time Taken")
    #
    # # This graph is for the total number of compilation warnings and it uses the same code as the above plot but with different data
    # p2 = bplt.figure(y_range=(0, max(compile_warn) + 50), plot_width=1000, plot_height=600,
    #                  tools="tap,wheel_zoom,box_zoom,reset",
    #                  title="Compilation Warnings", y_axis_label="Number of Compilation Warnings", x_axis_label="Date",
    #                  x_axis_type="datetime",
    #                  toolbar_location="below", sizing_mode='scale_width')
    # p2.circle('nt', 'cmt', size=10, color="blue", source=src)
    # p2.line('nt', 'cmt', color="blue", source=src)
    # taptool = p2.select(type=btools.TapTool)
    # taptool.callback = bcall.OpenURL(url=url)
    # tab3 = Panel(child=p2, title="Compile Time Warnings")
    #
    # src1 = bplt.ColumnDataSource(data=dict(cts=counts,
    #                                        wts=warning_types_list))
    #
    # # This plot is a bar graph showing the top 7 thorns with the most warnings
    # p3 = bplt.figure(x_range=warning_types_list, plot_width=1200, title="Compilation Warning Thorns",
    #                  y_axis_label="Number of Warnings", x_axis_label="Name of Thorn",
    #                  toolbar_location="below", tools="tap,wheel_zoom,box_zoom,reset", sizing_mode='scale_width')
    # p3.vbar(x='wts', top='cts', width=0.9, source=src1,
    #         line_color='white', fill_color=factor_cmap('wts', palette=viridis(len(counts)), factors=warning_types_list))
    # tab4 = Panel(child=p3, title="Compilation Warning Thorns")
    #
    # p.xaxis.major_label_orientation = math.pi / 6
    # p1.xaxis.major_label_orientation = math.pi / 6
    # p2.xaxis.major_label_orientation = math.pi / 6

    # Bokeh createst the html script and javscript for the plots using this code
    html = file_html(Tabs(tabs=[tab1]), CDN, "Plots")
    with open("./docs/plot.html", "w") as fp:
        fp.write(html)
    #script, div = components(Tabs(tabs=[tab1]))
    script, div = components(p)
    return script, div



def summary_to_html(readfile, writefile):
    '''
        This function reads the log file and outputs and html
        page with the summary in a table
    '''

    data = create_summary(readfile)
    print(data)
    contents = ""
    script, div = plot_test_data(readfile)

    # Check Status Using the data from the summary
    status = "All Tests Passed"
    if data["Number of tests passed"] == 0:
        status = "All Tests Passed"
    elif data["Number failed"] != 0:
        status = "Some Tests Failed"
        # Send email if tests failed
        # os.system(f'python3 mail.py')
    with open(writefile, "w") as fp:
        for key in data.keys():
            # Add a table row for each data field
            contents += f"        <tr><th>{key}</th><td>{data[key]}</td><tr>\n"

        # The formatted string holds the html template and loads in the values for content and status
        template = f'''<!doctype html>
    <html lang="en">
        <head>
            <title>Summary of Tests</title>
            <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css" integrity="sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T" crossorigin="anonymous">
            <style>
            .bk-root .bk {{
                margin: 0 auto !important;
            }}
            </style>
            <style>
            .sidebar {{
                height: 100%; 
                width: 150px;
                position: fixed;
                z-index: 1; 
                top: 0; 
                left: 0;
                background-color: #212529; 
                overflow-x: hidden;
                padding-top: 20px; 
            }}
            .sidebar a {{
                padding: 6px 8px 6px 16px;
                text-decoration: none;
                font-size: 18px;
                color: #dbdcdd;
                display: block;
                }}
            .sidebar a:hover {{
                color: white;
            }}
            .container{{
              padding-left: 150px;
              font-size: 18px;
            }}
                        /* On screens that are less than 700px wide, make the sidebar into a topbar */
            @media screen and (max-width: 500px) {{
            .sidebar {{
              display: none;
            }}
            .container {{
              padding-left:0px;
            }}
            }}
            </style>
            <script src="https://cdn.bokeh.org/bokeh/release/bokeh-2.0.1.min.js"
            crossorigin="anonymous"></script>
            {script}

        </head>
        <body>
            <div class="sidebar">
            </div>
            <script src='version.js'>
            </script>
            <div class="container">
                <h1 style="text-align:center">{status}</h1>
                <table class="table table-bordered " >
                <caption style="text-align:center;font-weight: bold;caption-side:top">Summary</caption>
                {contents}
                <br>
                <table class="table table-bordered " >
                <caption style="text-align:center;font-weight: bold;caption-side:top">Commits in Last Push</caption>
                {gen_commits()}
                </table>
                {gen_diffs(readfile)}
                <br>
                <br>
                <table style="margin: 0 auto;">
                    <iframe src="plot.html" height="700" width="1100"></iframe>
                </table>
                <table style="margin: 0 auto;">
                    {div}
                </table>
            <div>

        </body>
    </html>
        '''
        fp.write(template)


def write_to_csv(readfile):
    '''
        This function is used to store data between builds into a csv
    '''

    # total = sum(x[1] for x in get_times(readfile).items())  #

    data = create_summary(readfile)
    print("hello")
    #data["Time Taken"] = total / 60
    local_time = str(int(time.time()))  # datetime.today().strftime('%s') to convert to unix timestamp instead of a
    # normal date format. This helps in plotting as all x axis elements are now unique
    # local_time+=f"({curr_ver})"
    #data["Compile Time Warnings"] = get_compile(f"records/version_{curr_ver}/build_{curr_ver}.log")
    with open('new_test_nums.csv', 'a') as csvfile:
        contents = f"{local_time}"
        for key in data.keys():
            contents += f",{data[key]}"
        contents += f",{curr_ver}"
        contents += "\n"
        csvfile.write(contents)


# import glob, glob.glob("records/*/"build_1_2_*")


if __name__ == "__main__":
    write_to_csv(curr)
    summary_to_html(curr, "docs/new_index.html")
    copy_index(get_version() - 1)
    test_comparison = test_comp(curr, last)
    if len(test_comparison["Failed Tests"]) != 0 or len(test_comparison["Newly Passing Tests"]) != 0:
        os.system("python3 mail.py")