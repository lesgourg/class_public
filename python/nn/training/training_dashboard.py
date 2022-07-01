import curses
from tabulate import tabulate

def list_of_dicts_to_dict_of_lists(l):
    if not l:
        return {}
    keys = l[0].keys()
    return {k: [d[k] for d in l] for k in keys}

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

class Dashboard:
    pass

class SimpleDashboard(Dashboard):
    def __init__(self):
        self.table_style = "psql"

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def update(self, phase, epoch, epochs, batch, batches, loading_time, live, history):
    # def print_info(self, phase, epoch, i, batches, running_loss, timer):
        print("TRAINING PHASE:", phase)
        print("EPOCH: {:3d} / {:3d}".format(epoch, epochs))
        print("BATCH: {:5d} / {:5d}".format(batch, batches))
        print()

        print("LIVE METRICS:")
        print("=============")
        print("* Data Loading / cosmology: {:.3f} ms".format(loading_time * 1e3))

        live = list_of_dicts_to_dict_of_lists(live)
        live["loss"] = ["{:.03e}".format(v) for v in live["loss"]]
        live["forward"] = ["{:.03f} ms".format(v * 1e3) for v in live["forward"]]
        live["backward"] = ["{:.03f} ms".format(v * 1e3) for v in live["backward"]]

        table = tabulate(
            live, headers="keys",
            tablefmt=self.table_style, disable_numparse=True)
        print(table)
        print()
        print()

class CursesDashboard(Dashboard):

    def __init__(self):
        self.bar_width = 40
        self.table_style = "psql"
        self.reset_line()

    def setup(self):
        self.scr = curses.initscr()
        curses.cbreak()
        curses.noecho()
        self.scr.keypad(True)

    def teardown(self):
        curses.nocbreak()
        self.scr.keypad(False)
        curses.echo()
        curses.endwin()

    def __enter__(self):
        self.setup()
        return self

    def __exit__(self, *args):
        self.teardown()

    def reset_line(self):
        self.current_line = 0

    def add_line(self, line, column=0):
        self.scr.addstr(self.current_line, column, line)
        self.current_line += 1 + line.count("\n")

    def add_line_underlined(self, line, char="-"):
        self.add_line(line)
        self.add_line(char*len(line))

    def skip_line(self):
        self.current_line += 1

    def update(self, phase, epoch, epochs, batch, batches, loading_time, live, history):
        self.scr.clear()

        # 0) Header
        self.reset_line()
        phase_string = "Phase: {}".format(phase)
        self.add_line_underlined(phase_string)

        # 1) Progress
        self.add_line(self.get_epoch_bar(epoch, epochs))
        self.add_line(self.get_batch_bar(batch, batches))
        self.skip_line()

        # 2) Live Info
        self.update_live_info(loading_time, live)
        self.skip_line()

        # 3) History
        self.update_history(history)

        self.scr.refresh()

    def update_live_info(self, loading_time, live):
        self.add_line_underlined("LIVE METRICS:")
        self.add_line("* Data Loading / cosmology: {:.3f} ms".format(loading_time * 1e3))
        self.skip_line()

        live = list_of_dicts_to_dict_of_lists(live)
        live["loss"] = ["{:.03e}".format(v) for v in live["loss"]]
        live["forward"] = ["{:.03f} ms".format(v * 1e3) for v in live["forward"]]
        live["backward"] = ["{:.03f} ms".format(v * 1e3) for v in live["backward"]]

        table = tabulate(
            live, headers="keys",
            tablefmt=self.table_style, disable_numparse=True)
        row_count = table.count("\n") + 1
        self.add_line(table)

    def update_history(self, history, limit=4):
        if len(history) > limit:
            keys = history.keys()
            show_title = True
            for key_chunk in chunks(sorted(keys), limit):
                subhist = {k: history[k] for k in key_chunk}
                self.update_history_single(subhist, show_title=show_title)
                self.skip_line()
                show_title = False
        else:
            self.update_history_single(history)


    def update_history_single(self, history, show_title=True):
        if show_title:
            self.add_line_underlined("VALIDATION LOSS HISTORY:")

        history = {k: ["{:.03e}".format(v) for v in history[k]] for k in history}
        if history:
            index = list(range(1, 1 + len(history[next(iter(history.keys()))])))
        else:
            index = []

        history_table = tabulate(
            history, headers="keys", tablefmt=self.table_style,
            disable_numparse=True, showindex=index)
        self.add_line(history_table)

    def create_progress_bar(self, progress, width):
        assert width > 2
        avail = width - 2
        filled = int(progress * avail)
        empty = avail - filled
        arrow_head = ">" if filled >= 1 else ""
        return "[" + ("=" * (max(0, filled - 1))) + arrow_head + (" " * empty) + "]"

    def get_epoch_bar(self, epoch, epochs):
        bar = self.create_progress_bar(epoch / epochs, width=self.bar_width)
        return bar + " EPOCH {} / {}".format(epoch + 1, epochs)

    def get_batch_bar(self, batch, batches):
        bar = self.create_progress_bar(batch / batches, width=self.bar_width)
        return bar + " BATCH {} / {}".format(batch + 1, batches)

