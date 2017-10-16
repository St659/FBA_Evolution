import unittest
import xlrd
import csv
from FBA_Evolver import FBAEvolver, FBAPlotter
from matplotlib import pyplot as plt
import numpy as np
import os

class TestFBAEvolver(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._evolver = FBAEvolver('E:\\Chrome Download\\FBA\\FBA\\iPRAE34_edit.xls')


    def test_FBA_active(self):
        self.assertEqual(self._evolver.reader_active, True,)

    def test_num_exchange(self):
        #num_reactions, exchange_start, o2_pos = self._evolver.find_number_of_exchanges()
        #self.assertEqual(num_reactions, 331)
        #self.assertEqual(exchange_start + num_reactions -1,709 )

        self.assertEqual(self._evolver.start, 379)
        self.assertEqual(self._evolver.o2, 628)

    def test_set_exchange_values(self):
        set_values = [0]*331
        set_values[0] = 1
        oxygen = 628-379
        set_values[oxygen] = 1
        print(set_values)

        self._evolver.set_exchanges(set_values)
        workbook = xlrd.open_workbook('E:\\Chrome Download\\FBA\\FBA\\iPRAE34_edit.xls')
        ws = workbook.sheet_by_name('Reactions')
        nrows = ws.nrows
        rows = ws.col_slice(2, start_rowx=378, end_rowx=709)
        names = ws.col_slice(0, start_rowx=378, end_rowx=709)
        self.assertEqual(names[0].value, 'EX_12ppd_R(e)')
        self.assertEqual(rows[0].value, -1000)
        self.assertEqual(names[-1].value, 'EX_zn2(e)')
        self.assertEqual(rows[-1].value, 0)
        #Test oxygen level is still - 20
        self.assertEqual(names[oxygen].value, 'EX_o2(e)')
        self.assertEqual(rows[oxygen].value, -20)




    def test_FBA_runs(self):
        growth = self._evolver.run_fba(self._evolver.model, 'E:\\Chrome Download\\FBA\\FBA\\out.xls')
        self.assertIsNotNone(growth)


class TestFBAPlotter(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        file_0 = [0]*331
        file_1 = [0]*331

        file_0[0] = 1
        file_0[1] = 1
        print(file_0)

        file_1[1] = 1
        cls._filenames = list()

        files = [0,1]
        for file, data in zip(files, [file_0, file_1]):
            cls._filename = '.\\Test_Results\\' + str(file) +'.csv'
            with open(cls._filename, 'w') as csvfile:
                writer = csv.writer(csvfile, delimiter=',')
                writer.writerow(data)
                cls._filenames.append(cls._filename)

    def test_read_files(self):
        filenames = list()
        for f_file in os.listdir('./26-5-17'):
            filenames.append('./26-5-17/' + f_file)
        plotter_one = FBAPlotter(filenames, 'E:\\Chrome Download\\FBA\\FBA\\iPRAE34_plot_test.xls')
        self.assertEqual(len(plotter_one.fitnesses), 5)

    def test_number_of_files(self):
        plotter_one = FBAPlotter(self._filename, 'E:\\Chrome Download\\FBA\\FBA\\iPRAE34_plot_test.xls')
        plotter_two = FBAPlotter(self._filenames, 'E:\\Chrome Download\\FBA\\FBA\\iPRAE34_plot_test.xls')

        self.assertEqual(plotter_one.num_files, 1)
        self.assertEqual(plotter_two.num_files, 2)

    def test_get_set_reactions(self):
        plotter = FBAPlotter(self._filenames, 'E:\\Chrome Download\\FBA\\FBA\\iPRAE34_plot_test.xls')
        self.assertEqual(len(plotter.unique), 2)
        self.assertEqual(plotter.unique[0], 'EX_12ppd_R(e)')

    def test_filter_reaction_count(self):
        plotter = FBAPlotter(self._filenames, 'E:\\Chrome Download\\FBA\\FBA\\iPRAE34_plot_test.xls')
        filtered_name, filtered_count = plotter.filter_reactions(1)
        self.assertEqual(len(filtered_name), 1)

    def test_number_reactions(self):
        plotter = FBAPlotter(self._filenames, 'E:\\Chrome Download\\FBA\\FBA\\iPRAE34_plot_test.xls')
        self.assertEqual(plotter.active_exchange_count[0], 1)
        self.assertEqual(plotter.active_exchange_count[1], 2)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ind = np.arange(len(plotter.unique))
        width = 0.35
        ax.bar(ind, plotter.active_exchange_count, width)
        ax.set_xticks(ind + width/ 2)
        ax.set_xticklabels(plotter.unique)
        plt.show()