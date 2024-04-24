import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
from scipy import signal, stats
import numpy as np 
import random
import pickle
import csv

COLOR_LIST = list(mcolors.XKCD_COLORS)

class AnalyzeNeurons():
    """Analyze segmented calcium recordings."""
    def __init__(self):
        pass

    def find_file(self, folder_path: str, target: str) -> list[str]:
        """Find all dff csv files in the given folder."""
        file_list = []
        for (folder_path, _, fnames) in os.walk(folder_path):
            for file_name in fnames:
                if file_name == target or file_name.endswith(target):
                    file_path = os.path.join(folder_path, file_name)
                    file_list.append(file_path)

        return file_list

    def read_csv(self, path: str) -> pd.DataFrame:
        """Read the csv file."""
        with open(path, "r") as file:
            dff_file = pd.read_csv(file)
        return dff_file

    def _analyze(self, roi_dict: dict | None, cell_size: dict, 
                 dff: dict | pd.DataFrame, mda_path: str, save_path: str):
        """Analyze Cells."""
        roi_dff, n_dff, spk = self._analyze_dff(dff) 
        framerate, binning, pixel_size, objective, total_frames = self._extract_metadata(mda_path)

        if not cell_size:
            cell_size, cs_arr = self._cell_size(roi_dict, binning, pixel_size, objective)
        
        roi_analysis = self._analyze_roi(roi_dff, spk, framerate)
        mean_connect = self._get_mean_connect(roi_dff, spk)

        self._save_results(save_path, spk, cell_size, roi_analysis, framerate, total_frames)
        self._generate_summary(save_path, roi_analysis, spk, '/summary.txt', cell_size,
                                framerate, total_frames, mean_connect)

        # plot calcium traces
        self._plot_traces(n_dff, spk, save_path)

    def _reanalyze(self, folder_path: str, save_path: str):
        """Renalyze the data."""
        # traverse through the folder selected, go to one folder
        dff_file = os.path.join(folder_path, "dff.csv")
        mda_file = folder_path + "_metadata.txt"

        if os.path.exists(dff_file) and os.path.exists(mda_file):
            # get cell size
            path_1 = os.path.join(folder_path, "roi_size.csv")
            path_2 = os.path.join(folder_path, "roi_data.csv")

            if os.path.exists(path_1):
                cell_size = self._get_cell_size(path_1, "cell size")
            elif os.path.exists(path_2):
                cell_size = self._get_cell_size(path_2, "cell_size (um)")    
        
            dff_df = self.read_csv(dff_file)
            self._analyze(None, cell_size, dff_df, mda_file, save_path)

        else:
            print("one of the required files (dff.csv or metadata.txt) not found. Please run segmentation with the analysis.")

    def _analyze_dff(self, dff_df: pd.DataFrame | dict) -> tuple[dict, dict, dict]:
        """Analyze the dff file and get spikes."""
        if isinstance(dff_df, dict):
            dff_col = dff_df.keys()
        else:
            dff_col = dff_df.columns

        roi_dff = {}
        spk_times = {}
        n_dff = {}
        for col in dff_col:
            dff = list(dff_df[col]) # roi_dff[col] is the dff of one ROI -> LIST
            roi_dff[int(col)]= dff
            nn_dff = (dff - np.min(dff)) / (np.max(dff) - np.min(dff))
            n_dff[int(col)] = nn_dff
            spikes = self._mean_spk(roi_dff[int(col)], 90, 0.25)
            spk_times[int(col)] = list(spikes)

        return roi_dff, n_dff, spk_times

    def _cell_size(self, roi_dict: dict, binning: int, pixel_size: int, 
                   objective: int) -> tuple[dict, np.ndarray]:
        """Calculate cell size."""
        cs_dict = {}
        for r in roi_dict:
            cs_dict[r] = len(roi_dict[r]) * binning * pixel_size / objective # pixel to um

        cs_arr = np.array(list(cs_dict.items()))

        return cs_dict, cs_arr

    def _mean_spk(self, roi_dff_list: list, frame_window_ptg: int, prom_pctg: float) -> dict[dict]:
        """Test threshold for the scipy find peaks for one ROI."""
        start_frame = len(roi_dff_list) - int(len(roi_dff_list)*(frame_window_ptg/100))
        # using different percentage of the median
        prominence = np.mean(roi_dff_list[start_frame:]) * prom_pctg
        peaks, _ = signal.find_peaks(roi_dff_list, prominence=prominence)
        spk_time_ptg = list(peaks)

        return spk_time_ptg

    def _extract_metadata(self, mda_file: str)-> tuple[float, int, float, int, int]:
        """Extract information from the metadata."""
        framerate = 0
        fr = False
        bn = False
        obj = False
        ps = False
        t_f = False

        with open(mda_file) as f:
            metadata = f.readlines()

        for line in metadata:
            line = line.strip()
            if line.startswith('"Exposure-ms": '):
                exposure = float(line[15:-1]) / 1000  # exposure in seconds
                framerate = 1 / exposure  # frames/second
                fr = True
            if line.startswith('"Binning": '):
                binning = int(line[11:-1])
                bn = True
            if line.startswith('"PixelSizeUm": '):
                pixel_size = float(line[15:-1])
                ps = True
            if line.startswith('"Nosepiece-Label": '):
                words = line[19:-1].strip('\"').split(" ")
                objective = int([word for word in words if word.endswith("x")][0][:-1])
                obj = True
            if line.startswith('"Frames": '):
                total_frames = int(line[len('"Frames": '):-1])
            if fr and bn and ps and obj and t_f:
                break

        return framerate, binning, pixel_size, objective, total_frames

    def _analyze_roi(self, roi_dff: dict, spk_times: dict, framerate: float):
        """Analyze the dff from each ROI."""
        amplitude_info = self._get_amplitude(roi_dff, spk_times)
        time_to_rise = self._get_time_to_rise(amplitude_info, framerate)
        max_slope = self._get_max_slope(roi_dff, amplitude_info)
        iei = self._analyze_iei(spk_times, framerate)
        roi_analysis = amplitude_info

        for r in roi_analysis:
            # roi_analysis[r]['spike_times'] = spk_times[r]
            roi_analysis[r]['time_to_rise'] = time_to_rise[r]
            roi_analysis[r]['max_slope'] = max_slope[r]
            roi_analysis[r]['IEI'] = iei[r]

        return roi_analysis

    def _get_amplitude(self, roi_dff: dict, spk_times: dict,
                        deriv_threhold=0.01, reset_num=17, neg_reset_num=2, total_dist=40) -> dict:
        """Get the amplitude."""
        amplitude_info = {}

        # for each ROI
        for r in spk_times:
            amplitude_info[r] = {}
            amplitude_info[r]['amplitudes'] = []
            amplitude_info[r]['peak_indices'] = []
            amplitude_info[r]['base_indices'] = []

            if len(spk_times[r]) > 0:
                dff_deriv = np.diff(roi_dff[r]) # the difference between each spike

                # for each spike in the ROI
                for i in range(len(spk_times[r])):
                    # Search for starting index for current spike
                    searching = True
                    under_thresh_count = 0
                    total_count = 0
                    start_index = spk_times[r][i] # the frame for the first spike

                    if start_index > 0:
                        while searching:
                            start_index -= 1
                            total_count += 1

                            # If collide with a new spike
                            if start_index in spk_times[r]:
                                subsearching = True
                                negative_count = 0

                                while subsearching:
                                    start_index += 1
                                    if start_index < len(dff_deriv):
                                        if dff_deriv[start_index] < 0:
                                            negative_count += 1

                                        else:
                                            negative_count = 0

                                        if negative_count == neg_reset_num:
                                            subsearching = False
                                    else:
                                        subsearching = False

                                break

                            # if the difference is below threshold
                            if dff_deriv[start_index] < deriv_threhold:
                                under_thresh_count += 1
                            else:
                                under_thresh_count = 0

                            # stop searching for starting index
                            if under_thresh_count >= reset_num or start_index == 0 or total_count == total_dist:
                                searching = False

                    # Search for ending index for current spike
                    searching = True
                    under_thresh_count = 0
                    total_count = 0
                    end_index = spk_times[r][i]

                    if end_index < (len(dff_deriv) - 1):
                        while searching:
                            end_index += 1
                            total_count += 1

                            # If collide with a new spike
                            if end_index in spk_times[r]:
                                subsearching = True
                                negative_count = 0
                                while subsearching:
                                    end_index -= 1
                                    if dff_deriv[end_index] < 0:
                                        negative_count += 1
                                    else:
                                        negative_count = 0
                                    if negative_count == neg_reset_num:
                                        subsearching = False
                                break
                            if dff_deriv[end_index] < deriv_threhold:
                                under_thresh_count += 1
                            else:
                                under_thresh_count = 0

                            # NOTE: changed the operator from == to >=
                            if under_thresh_count >= reset_num or end_index >= (len(dff_deriv) - 1) or \
                                    total_count == total_dist:
                                searching = False

                    # Save data
                    spk_to_end = roi_dff[r][spk_times[r][i]:(end_index + 1)]
                    start_to_spk = roi_dff[r][start_index:(spk_times[r][i] + 1)]
                    try:
                        amplitude_info[r]['amplitudes'].append(np.max(spk_to_end) - np.min(start_to_spk))
                        amplitude_info[r]['peak_indices'].append(int(spk_times[r][i] + np.argmax(spk_to_end)))
                        amplitude_info[r]['base_indices'].append(int(spk_times[r][i] -
                                                                    (len(start_to_spk) - (np.argmin(start_to_spk) + 1))))
                    except ValueError:
                        pass

        return amplitude_info

    def _get_time_to_rise(self, amplitude_info: dict, framerate: float) -> dict:
        """Get time to rise"""
        time_to_rise = {}
        for r in amplitude_info:
            time_to_rise[r] = []
            if len(amplitude_info[r]['peak_indices']) > 0:
                for i in range(len(amplitude_info[r]['peak_indices'])):
                    peak_index = amplitude_info[r]['peak_indices'][i]
                    base_index = amplitude_info[r]['base_indices'][i]
                    frames = peak_index - base_index + 1
                    if framerate:
                        time = frames / framerate  # frames * (seconds/frames) = seconds
                        time_to_rise[r].append(time)
                    else:
                        time_to_rise[r].append(frames)

        return time_to_rise

    def _get_max_slope(self, roi_dff: dict, amplitude_info: dict):
        """Get Max slope"""
        max_slope = {}
        for r in amplitude_info:
            max_slope[r] = []
            dff_deriv = np.diff(roi_dff[r])
            if len(amplitude_info[r]['peak_indices']) > 0:
                for i in range(len(amplitude_info[r]['peak_indices'])):
                    peak_index = amplitude_info[r]['peak_indices'][i]
                    base_index = amplitude_info[r]['base_indices'][i]
                    slope_window = dff_deriv[base_index:(peak_index + 1)]
                    max_slope[r].append(np.max(slope_window))

        return max_slope

    def _analyze_iei(self, spk_times: dict, framerate: float):
        """Analyze IEI."""
        iei = {}
        for r in spk_times:
            iei[r] = []

            if len(spk_times[r]) > 1:
                iei_frames = np.diff(np.array(spk_times[r]))
                if framerate:
                    iei[r] = iei_frames / framerate # in seconds
                else:
                    iei[r] = iei_frames
        return iei

    def _get_mean_connect(self, roi_dff: dict, spk_times: dict):
        """Calculate connectivity."""
        A = self._get_connect_matrix(roi_dff, spk_times)

        if A is not None:
            if len(A) > 1:
                mean_connect = np.median(np.sum(A, axis=0) - 1) / (len(A) - 1)
            else:
                mean_connect = 'N/A - Only one active ROI'
        else:
            mean_connect = 'No calcium events detected'

        return mean_connect
    
    def _get_connect_matrix(self, roi_dff: dict, spk_times: dict) -> np.ndarray:
        """Calculate the connectivity matrix."""
        active_roi = [r for r in spk_times if len(spk_times[r]) > 0]

        if len(active_roi) > 0:
            phases = {}
            for r in active_roi:
                phases[r] = self._get_phase(len(roi_dff[r]), spk_times[r])

            connect_matrix = np.zeros((len(active_roi), len(active_roi)))
            for i, r1 in enumerate(active_roi):
                for j, r2 in enumerate(active_roi):
                    connect_matrix[i, j] = self._get_sync_index(phases[r1], phases[r2])
        else:
            connect_matrix = None

        return connect_matrix
    
    def _get_phase(self, total_frames: int, spks: list) -> list:
        """Get Phase"""
        spikes = spks.copy()
        if len(spikes) == 0 or spikes[0] != 0:
            spikes.insert(0, 0)
        if spikes[-1] != (total_frames - 1):
            spikes.append(total_frames - 1)

        phase = []
        for k in range(len(spikes) - 1):
            t = spikes[k]

            while t < spikes[k + 1]:
                instant_phase = (2 * np.pi) * ((t - spikes[k]) / \
                                               (spikes[k+1] - spikes[k])) + (2 * np.pi * k)
                phase.append(instant_phase)
                t += 1
        phase.append(2 * np.pi * (len(spikes) - 1))

        return phase

    def _get_sync_index(self, x_phase: list, y_phase: list):
        """Calculate the pair-wise synchronization index of the two ROIs"""
        phase_diff = self._get_phase_diff(x_phase, y_phase)
        sync_index = np.sqrt((np.mean(np.cos(phase_diff)) ** 2) + (np.mean(np.sin(phase_diff)) ** 2))

        return sync_index
    
    def _get_phase_diff(self, x_phase: list, y_phase: list) -> np.ndarray:
        """Calculate the absolute phase difference between two calcium
            traces x and y phases from two different ROIs."""
        x_phase = np.array(x_phase)
        y_phase = np.array(y_phase)
        phase_diff = np.mod(np.abs(x_phase - y_phase), (2 * np.pi))

        return phase_diff

    def _plot_traces(self, roi_dff: dict, spk_times: dict, save_path: str) -> None:
        """Plot  traces."""
        dff_to_plot, color_list = self._random_pick(roi_dff, 10)
        self._plot_traces_no_peaks(roi_dff, dff_to_plot, color_list, save_path)
        self._plot_traces_w_peaks(roi_dff, dff_to_plot, spk_times, color_list, save_path)

    def _plot_traces_no_peaks(self, roi_dff: dict, dff_to_plot: list,
                              color_list: list, path: str) -> None:
        """Plot the traces."""
        fig, ax = plt.subplots(figsize=(20, 20))
        if len(dff_to_plot) > 0:
            dff_max = np.zeros(len(dff_to_plot))
            for max_index, dff_index in enumerate(dff_to_plot):
                dff_max[max_index] = np.max(roi_dff[dff_index])
            height_increment = max(dff_max)

            y_pos = []
            for height_index, d in enumerate(dff_to_plot):
                y_pos.append(height_index * (1.2 * height_increment))
                ax.plot(roi_dff[d] + height_index * (1.2 * height_increment), color=color_list[height_index],
                        linewidth=3)
            ax.set_yticks(y_pos, labels=dff_to_plot)
            fname = "traces_no_detection.png"
            plt.savefig(os.path.join(path, fname))
            plt.close()

    def _plot_traces_w_peaks(self, roi_dff: dict, dff_to_plot: list, spk_times: dict, 
                             color_list: list, path: str):
        """Plot traces with peak detected."""
        fig, ax = plt.subplots(figsize=(20, 20))
        if len(dff_to_plot) > 0:
            dff_max = np.zeros(len(dff_to_plot))
            for max_index, dff_index in enumerate(dff_to_plot):
                dff_max[max_index] = np.max(roi_dff[dff_index])
            height_increment = max(dff_max)

            y_pos = []
            for height_index, d in enumerate(dff_to_plot):
                y_pos.append(height_index * (1.2 * height_increment))
                ax.plot(roi_dff[d] + height_index * (1.2 * height_increment), color=color_list[height_index], linewidth=3)
                if len(spk_times[d]) > 0:
                    y = [roi_dff[d][i] for i in spk_times[d]]
                    ax.plot(spk_times[d], y + height_index * (1.2 * height_increment),
                                   ms=6, color='r', marker='o', ls='', label=f"{d}: {len(spk_times[d])}")
                    ax.legend()

            ax.set_yticks(y_pos, labels=dff_to_plot)
            fname = "traces_w_peaks.png"
            plt.savefig(os.path.join(path, fname))
            plt.close()

    def _random_pick(self, roi_dff: dict, num: int) -> tuple[list, list]:
        """Pick 10 traces randomly to plot."""
        num_f = np.min([num, len(roi_dff)])
        final_dff = random.sample(list(roi_dff.keys()), num_f)
        final_dff.sort()
        rand_color_ind = random.sample(COLOR_LIST, k=num_f)

        return final_dff, rand_color_ind
    
    def _save_results(self, save_path: str, spk: dict, cell_size: dict, roi_analysis: dict,
                      framerate: float, total_frames: int):
        """Save the analysis results."""
        # save spike times
        if not os.path.isdir(save_path):
            os.mkdir(save_path)

        with open(save_path + '/spike_times.pkl', 'wb') as spike_file:
            pickle.dump(spk, spike_file)

        # save 
        roi_data = self.all_roi_data(roi_analysis, cell_size, spk, framerate, total_frames)
        with open(save_path + '/roi_data.csv', 'w', newline='') as roi_data_file:
            writer = csv.writer(roi_data_file, dialect='excel')
            fields = ['ROI', 'cell_size (um)', '# of events', 'frequency (num of events/s)',
                    'average amplitude', 'amplitude SEM', 'average time to rise', 'time to rise SEM',
                    'average max slope', 'max slope SEM',  'InterEvent Interval', 'IEI SEM']
            writer.writerow(fields)
            writer.writerows(roi_data)
            
    def all_roi_data(self, roi_analysis: dict, cell_size: dict, spk_times: dict, framerate: float, total_frames: int) -> list:
        """Compile data from roi_analysis into one csv file"""
        num_roi = len(roi_analysis.keys())
        if len(cell_size.keys()) == num_roi and len(spk_times.keys()) == num_roi:
            roi_data = np.zeros((num_roi, 12))
            recording_time = total_frames/framerate

            for i, r in enumerate(roi_analysis):
                roi_data[i, 0] = r
                roi_data[i, 1] = cell_size[int(r)]
                num_e = len(spk_times[r])
                roi_data[i, 2] = num_e
                roi_data[i, 3] = num_e / recording_time
                roi_data[i, 4] = np.mean(roi_analysis[r]['amplitudes'])
                roi_data[i, 5] = stats.sem(roi_analysis[r]['amplitudes'])
                # roi_data[i, 5] = roi_analysis[r]['amplitudes']
                roi_data[i, 6] = np.mean(roi_analysis[r]['time_to_rise'])
                roi_data[i, 7] = stats.sem(roi_analysis[r]['time_to_rise'])
                # roi_data[i, 7] = roi_analysis[r]['time_to_rise']
                roi_data[i, 8] = np.mean(roi_analysis[r]['max_slope'])
                roi_data[i, 9] = stats.sem(roi_analysis[r]['max_slope'])
                # roi_data[i, 9] = roi_analysis[r]['max_slope']
                roi_data[i, 10] = np.mean(roi_analysis[r]['IEI'])
                roi_data[i, 11] = stats.sem(roi_analysis[r]['IEI'], nan_policy='omit')
        else:
            print('please make sure that the number of ROIs in each dictionary is the same')

        return roi_data

    def _get_cell_size(self, path: str, col: str) -> dict:
        """Extract cell size from previous analysis."""
        cell_data = self.read_csv(path)
        cell_size = {}
        rois = [roi for roi in cell_data["ROI"]]
        cs = [size for size in cell_data[col]]
            
        if len(rois) == len(cs):
            for i, roi in enumerate(rois):
                cell_size[roi] = cs[i]

        return cell_size

    def _generate_summary(self, save_path:str, roi_analysis: dict,
                         spike_times: dict, file_name: str,
                         cell_size: dict, framerate: float, total_frames: int, 
                         mean_connect: dict) -> None:
        """Generate summary."""
        avg_cs = float(np.mean(list(cell_size.values())))
        std_cs = float(np.std(list(cell_size.values())))

        total_amplitude = []
        total_time_to_rise = []
        total_max_slope = []
        total_IEI = []
        total_num_events = []

        for r in roi_analysis:
            if len(roi_analysis[r]['amplitudes']) > 0:
                total_amplitude.extend(roi_analysis[r]['amplitudes'])
                total_time_to_rise.extend(roi_analysis[r]['time_to_rise'])
                total_max_slope.extend(roi_analysis[r]['max_slope'])
                if len(spike_times[r]) > 0:
                    total_num_events.append(len(spike_times[r]))
            if len(roi_analysis[r]['IEI']) > 0:
                total_IEI.extend(roi_analysis[r]['IEI'])

        if any(spike_times.values()):
            avg_amplitude = np.mean(np.array(total_amplitude))
            std_amplitude = np.std(np.array(total_amplitude))
            avg_max_slope = np.mean(np.array(total_max_slope))
            std_max_slope = np.std(np.array(total_max_slope))
            # get the average num of events
            avg_num_events = np.mean(np.array(total_num_events))
            std_num_events = np.std(np.array(total_num_events))
            avg_time_to_rise = np.mean(np.array(total_time_to_rise))
            # avg_time_to_rise = f'{avg_time_to_rise} {units}'
            std_time_to_rise = np.std(np.array(total_time_to_rise))
            if len(total_IEI) > 0:
                avg_IEI = np.mean(np.array(total_IEI))
                # avg_IEI = f'{avg_IEI} {units}'
                std_IEI = np.std(np.array(total_IEI))
            else:
                avg_IEI = 'N/A - Only one event per ROI'
        else:
            avg_amplitude = 'No calcium events detected'
            avg_max_slope = 'No calcium events detected'
            avg_time_to_rise = 'No calcium events detected'
            avg_IEI = 'No calcium events detected'
            avg_num_events = 'No calcium events detected'
        percent_active = self._analyze_active(spike_times)

        with open(save_path + file_name, 'w') as sum_file:
            units = "seconds" if framerate else "frames"
            _, filename = os.path.split(save_path)
            sum_file.write(f'File: {filename}\n')
            if framerate:
                sum_file.write(f'Framerate: {framerate} fps\n')
            else:
                sum_file.write('No framerate detected\n')
            sum_file.write(f'Total ROI: {len(cell_size.keys())}\n')
            sum_file.write(f'Percent Active ROI (%): {percent_active}\n')

            # NOTE: include cell size in the summary text file
            sum_file.write(f'Average Cell Size(um): {avg_cs}\n')
            sum_file.write(f'\tCell Size Standard Deviation: {std_cs}\n')

            sum_file.write(f'Average Amplitude: {avg_amplitude}\n')

            if len(total_amplitude) > 0:
                sum_file.write(f'\tAmplitude Standard Deviation: {std_amplitude}\n')
            sum_file.write(f'Average Max Slope: {avg_max_slope}\n')
            if len(total_max_slope) > 0:
                sum_file.write(f'\tMax Slope Standard Deviation: {std_max_slope}\n')
            sum_file.write(f'Average Time to Rise ({units}): {avg_time_to_rise}\n')
            if len(total_time_to_rise) > 0:
                sum_file.write(f'\tTime to Rise Standard Deviation: {std_time_to_rise}\n')
            sum_file.write(f'Average Interevent Interval (IEI) ({units}): {avg_IEI}\n')
            if len(total_IEI) > 0:
                sum_file.write(f'\tIEI Standard Deviation: {std_IEI}\n')
            sum_file.write(f'Average Number of events: {avg_num_events}\n')
            if len(total_num_events) > 0:
                sum_file.write(f'\tNumber of events Standard Deviation: {std_num_events}\n')
                if framerate:
                    frequency = avg_num_events/(total_frames/framerate)
                    sum_file.write(f'Average Frequency (num_events/s): {frequency}\n')

            sum_file.write(f'Global Connectivity: {mean_connect}')

    def _analyze_active(self, spk_times: dict):
        """Calculate the percentage of active cell bodies"""
        active = 0
        for r in spk_times:
            if len(spk_times[r]) > 0:
                active += 1
        active = (active / len(spk_times)) * 100
        return active