import tkinter
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
import pdbfinder as pf
import platform
import json
import os
import threading


OS = platform.system()
if OS == 'Darwin':
	separator = '/'
elif OS =='Windows':
	separator = '\\'
else:
	pass


mainwindow = tkinter.Tk()
mainwindow.title("Uniprot-PDB Mapping")
mainwindow.lift()
mainwindow.attributes('-topmost',True)
mainwindow.after_idle(mainwindow.attributes,'-topmost',False)

#Have the window open in the middle of the screen
screen_width = mainwindow.winfo_screenwidth()
screen_height = mainwindow.winfo_screenheight()

window_width = mainwindow.winfo_reqwidth()
window_height = mainwindow.winfo_reqheight()

positionRight = int(screen_width/2 - window_width/2)
positionDown = int(screen_height/2-window_height/2)

mainwindow.geometry("+{}+{}".format(positionRight, positionDown))


#Create the top half of the window containing the labels, text boxes and buttons
mainframe=tkinter.Frame(mainwindow)
mainframe.grid(row=0, column=0, sticky='NSEW')

topframe = tkinter.Frame(mainframe)
topframe.grid(row=0, column=0, sticky='NSEW', padx=5, pady=5)
sep = ttk.Separator(topframe, orient='horizontal')
sep.grid(row=1, column=0, sticky='NSEW')

#Create the labels
label1 = tkinter.Label(topframe, text="Input file:", anchor='w')
label1.grid(row=0, column=0, sticky=("NSEW"))

label2 = tkinter.Label(topframe, text="Output directory:", anchor='w')
label2.grid(row=1, column=0, sticky=('NSEW'))

label3 = tkinter.Label(topframe, text="Enter output file name:", anchor='w')
label3.grid(row=2, column=0, sticky=('NSEW'))

#Create the text boxes
ent1=tkinter.Entry(topframe, font=40)
ent1.grid(row=0, column=1, sticky=('NSEW'), padx=2, pady=2)

ent2 = tkinter.Entry(topframe, font=40)
ent2.grid(row=1, column=1, sticky=('NSEW'), padx=2, pady=2)

ent3 = tkinter.Entry(topframe, font=40)
ent3.grid(row=2, column=1, sticky=('NSEW'), padx=2, pady=2)


#Functions to control the buttons and the progress bar
def browsefile():
	filename = filedialog.askopenfilename(filetypes=(("Excel files", "*.xlsx"), ("All files", "*.*")))
	ent1.delete(0, tkinter.END)
	ent1.insert(tkinter.END, filename)

def browsedir(entry):
	dirname = filedialog.askdirectory()
	entry.delete(0, tkinter.END)
	entry.insert(tkinter.END, dirname)

def bar(pos):
	progbar['value'] = pos
	progbar.update()
	progbar.after(1)
	#print("Progbar update: %s" % pos)


def runprog(defaults):
	if os.path.isfile(ent1.get()) and os.path.exists(ent2.get()):
		if ent3.get() != '':
			threading.Thread(target=pf.process(ent1.get(), ent2.get() + separator, ent3.get(), defaults, bar, status)).start()
	#print(defaults)
		else:
			messagebox.showerror(title = "Error!", message = 'Output file name is blank.')
	else:
		messagebox.showerror(title = "Error!", message='Input file or output directory do not exist.')


def status(message = ""):
	status_label['text'] = message


def get_settings():
	print("Settings have been re-written")
	try:
		with open("settings.dat") as defaults_path:
			default_settings = json.load(defaults_path)
		return default_settings
	except:

		current_dir = os.getcwd()
		default_settings = {"eval_cutoff" : 20, "max_pdb" : 5, "min_percent" : 0, "query_length" : 200,
			"distance_cutoff" : 10, "pdb_path" : current_dir}

		return default_settings



default_settings = get_settings()


'''
Options/settings window
'''
def options_window(firstwindow, settings):



	def save_settings_to_var(settings):
		try:
			if int(float(blastq_entry.get())) < 1:
				raise ValueError
			else:
				settings["query_length"] = int(blastq_entry.get())
		except ValueError: #Thrown above and if the input is not a number
			messagebox.showerror(title="Error!", message="\"Blastp query length\" must be a positive integer.")
			blastq_entry.delete(0, tkinter.END)
			blastq_entry.insert(tkinter.END, settings["query_length"])
			return None

		try:
			if int(float(eval_entry.get())) < 1:
				raise ValueError
			else:
				settings["eval_cutoff"] = int(eval_entry.get())
		except ValueError:#Thrown above and if the input is not a number
			messagebox.showerror(title="Error!", message="\"Max e-value\" must be a positive integer.")
			eval_entry.delete(0, tkinter.END)
			eval_entry.insert(tkinter.END, settings['eval_cutoff'])
			return None

		try:
			if int(float(maxpdb_entry.get())) < 1:
				raise ValueError
			else:
				settings["max_pdb"] = int(maxpdb_entry.get())
		except ValueError:#Thrown above and if the input is not a number
			messagebox.showerror(title="Error!", message="\"Max number of PDB files to report\" must be a positive integer.")
			maxpdb_entry.delete(0, tkinter.END)
			maxpdb_entry.insert(tkinter.END, settings['max_pdb'])
			return None

		try:
			if int(float(minpercent_entry.get())) < 0:
				raise ValueError
			else:
				settings["min_percent"] = int(minpercent_entry.get())
		except ValueError:#Thrown above and if the input is not a number
			messagebox.showerror(title="Error!", message="\"Min percent identity\" must be ≥ 0.")
			minpercent_entry.delete(0, tkinter.END)
			minpercent_entry.insert(tkinter.END, settings['min_percent'])
			return None

		try:
			if float(distance_cutoff_entry.get()) <= 0.0:
				#print(distance_cutoff_entry.get())
				raise ValueError
			else:
				settings["distance_cutoff"] = float(distance_cutoff_entry.get())
		except ValueError as ve:#Thrown above and if the input is not a number
			print(ve)
			messagebox.showerror(title="Error!", message="\"Distance cutoff\" must be greater than 0.")
			distance_cutoff_entry.delete(0, tkinter.END)
			distance_cutoff_entry.insert(tkinter.END, settings["distance_cutoff"])
			return None

		if os.path.isdir(pdb_directory_entry.get()):
			settings["pdb_path"] = pdb_directory_entry.get()
		else:
			messagebox.showerror(title="Error!", message="\"PDB directory\" must specify a folder to store PDB files.")
			pdb_directory_entry.delete(0, tkinter.END)
			pdb_directory_entry.insert(tkinter.END, settings["pdb_path"])
			return None

		return True


	def apply_settings(frame, settings):
		if save_settings_to_var(settings):
			frame.grid_remove()
			mainwindow.title("Uniprot-PDB Mapping")

	def close_window():
		options.grid_remove()
		mainwindow.title("Uniprot-PDB Mapping")

	def save_settings(store):
		output_filename = filedialog.asksaveasfilename(title= "Select file", defaultextension ='.dat')
		save_settings_to_var(settings)
		with open(output_filename, 'w') as output:
			json.dump(store, output)

	def save_defaults():
		output_settings = {"query_length" : blastq_entry.get(), "eval_cutoff" : eval_entry.get(), "max_pdb" : maxpdb_entry.get(),
			"distance_cutoff" : distance_cutoff_entry.get(), "pdb_path" : pdb_directory_entry.get(), "min_percent" : minpercent_entry.get()}

		with open("settings.dat", 'w') as outfile:
			json.dump(output_settings, outfile)


	def load_settings(settings):
		settings_path = filedialog.askopenfilename(title="Select file", filetypes = (("settings files", "*.dat"), ("all files", "*.*")))

		with open(settings_path, 'r') as inputfile:
			loaded_settings = json.load(inputfile)

		try:
			for each in loaded_settings.keys():
				settings[each] = loaded_settings[each]

			blastq_entry.delete(0, tkinter.END)
			blastq_entry.insert(tkinter.END, settings["query_length"])
			eval_entry.delete(0, tkinter.END)
			eval_entry.insert(tkinter.END, settings['eval_cutoff'])
			maxpdb_entry.delete(0, tkinter.END)
			maxpdb_entry.insert(tkinter.END, settings['max_pdb'])
			minpercent_entry.delete(0, tkinter.END)
			minpercent_entry.insert(tkinter.END, settings['min_percent'])
			distance_cutoff_entry.delete(0, tkinter.END)
			distance_cutoff_entry.insert(tkinter.END, settings["distance_cutoff"])
			pdb_directory_entry.delete(0, tkinter.END)
			pdb_directory_entry.insert(tkinter.END, settings["pdb_path"])

		except:
			messagebox.showerror(title="Error!", message="Invalid file selection")

	settings_placeholder = {}

	options = tkinter.Frame(mainwindow)
	options.grid(row=0, column=0, sticky='NSEW')
	mainwindow.title("Options")

	options_top = tkinter.Frame(options)
	options_top.grid(row=0, column=0, sticky='NSEW', padx=5, pady=5)

	top_title = tkinter.Label(options_top, text="PDB finder settings:", font='Arial 14 bold')
	top_title.grid(row=0, column=0, columnspan=3, sticky='news')

	blastq_label = tkinter.Label(options_top, text="Blastp query length:", anchor='w')
	blastq_label.grid(row=1, column=0, sticky='NSEW')

	blastq_entry = tkinter.Entry(options_top, font=10)
	blastq_entry.grid(row=1, column=1, sticky='NSEW')
	blastq_entry.insert(tkinter.END, settings["query_length"])
	blastq_entry.bind("<Return>", (lambda event : apply_settings(options, settings)))

	eval_label = tkinter.Label(options_top, text="Max e-value to report:", anchor='w')
	eval_label.grid(row=2, column=0, sticky='NSEW')

	eval_entry = tkinter.Entry(options_top, font=10)
	eval_entry.grid(row=2, column=1, sticky='NSEW')
	eval_entry.insert(tkinter.END, settings["eval_cutoff"])
	eval_entry.bind("<Return>", (lambda event : apply_settings(options, settings)))

	maxpdb_label = tkinter.Label(options_top, text="Max number of PDB files to report", anchor='w')
	maxpdb_label.grid(row=3, column=0, sticky='NSEW')

	maxpdb_entry = tkinter.Entry(options_top, font=10)
	maxpdb_entry.grid(row=3, column=1, sticky='NSEW')
	maxpdb_entry.insert(tkinter.END, settings["max_pdb"])
	maxpdb_entry.bind("<Return>", (lambda event : apply_settings(options, settings)))

	minpercent_label = tkinter.Label(options_top, text="Min percent identity", anchor='w')
	minpercent_label.grid(row=4, column=0, sticky='NSEW')

	minpercent_entry = tkinter.Entry(options_top, font=10)
	minpercent_entry.grid(row=4, column=1, sticky='NSEW')
	minpercent_entry.insert(tkinter.END, settings["min_percent"])
	minpercent_entry.bind("<Return>", (lambda event : apply_settings(options, settings)))

	separator = ttk.Separator(options)
	separator.grid(row=1, column=0, sticky='NSEW')

	options_middle =tkinter.Frame(options)
	options_middle.grid(row=2, column=0, sticky='NSEW', padx=5, pady=3)

	middle_title = tkinter.Label(options_middle, text="Ligand finder settings", font='Arial 14 bold')
	middle_title.grid(row=0, column=0, columnspan=3, sticky='news')

	distance_cutoff_label = tkinter.Label(options_middle, text="Max distance to ligand (Å):", anchor='w')
	distance_cutoff_label.grid(row=1, column=0, sticky='NSEW')

	distance_cutoff_entry = tkinter.Entry(options_middle, font=10)
	distance_cutoff_entry.grid(row=1, column=1, sticky='NSEW')
	distance_cutoff_entry.insert(tkinter.END, settings["distance_cutoff"])
	distance_cutoff_entry.bind("<Return>", (lambda event : apply_settings(options, settings)))

	pdb_directory_label = tkinter.Label(options_middle, text="PDB Directory:", anchor='w')
	pdb_directory_label.grid(row=2, column=0, sticky='NSEW')

	pdb_directory_entry = tkinter.Entry(options_middle, font=10)
	pdb_directory_entry.grid(row=2, column=1, sticky='NSEW')
	pdb_directory_entry.insert(tkinter.END, settings["pdb_path"])
	pdb_directory_entry.bind("<Return>", (lambda event : apply_settings(options, settings)))

	pdb_directory_button = tkinter.Button(options_middle, text='Set directory', command=lambda: browsedir(pdb_directory_entry))
	pdb_directory_button.grid(row=2, column=2, sticky='NSEW')


	options_bottom = tkinter.Frame(options)
	options_bottom.grid(row=3, column=0, sticky="NSEW", padx=5, pady=3)

	applyb = tkinter.Button(options_bottom, text="Apply", command=lambda:apply_settings(options, settings))
	applyb.grid(row=1, column=1, sticky="NSEW")

	load_settings_button = tkinter.Button(options_middle, text='Load settings', command=lambda:load_settings(settings))
	load_settings_button.grid(row=3, column=0, sticky='NSEW')

	save_settings_button = tkinter.Button(options_middle, text='Save settings', command=lambda:save_settings(settings))
	save_settings_button.grid(row=3, column=2, sticky='NSEW')

	set_as_default = tkinter.Button(options_middle, text='Set as default settings', command=lambda:save_defaults())
	set_as_default.grid(row=3, column=1, sticky='NSEW')

	cancelb = tkinter.Button(options_bottom, text="Cancel", command=lambda:close_window())
	cancelb.grid(row=1, column=0, sticky="NSEW")

	options.grid_columnconfigure(0, weight=1)
	options.grid_rowconfigure(0, weight=1)

	options_top.grid_columnconfigure(0, weight=1)
	options_top.grid_columnconfigure(1, weight=10)


	options_top.grid_rowconfigure(0, weight=1)
	options_top.grid_rowconfigure(1, weight=1)
	options_top.grid_rowconfigure(2, weight=1)

	options_middle.grid_columnconfigure(0, weight=1)
	options_middle.grid_columnconfigure(1, weight=1)
	options_middle.grid_columnconfigure(2, weight=1)


	options_bottom.grid_columnconfigure(0, weight=1)
	options_bottom.grid_columnconfigure(1, weight=1)


#Create buttons
select_input = tkinter.Button(topframe, text="Select input", command=browsefile)
select_input.grid(row=0, column=2, sticky='NEWS')

select_output_dir = tkinter.Button(topframe, text="Select output directory", command=lambda:browsedir(ent2))
select_output_dir.grid(row=1, column=2, sticky='NEWS')

#Create the progress bar
progbar = ttk.Progressbar(mainframe, orient='horizontal', length=100, mode='determinate')
progbar.grid(row=3, column=0, sticky='NSEW', ipady=1)




#Create status label
status_label = tkinter.Label(mainframe, text='', anchor='w')
status_label.grid(row=2, column=0, sticky='NEWS')



#Create bottom button frame
bottomframe = tkinter.Frame(mainframe)
bottomframe.grid(row=4, column=0, sticky='NEWS', padx=5, pady=5)

#Create options button
options_button = tkinter.Button(bottomframe, text="Options", command=lambda:options_window(mainwindow, default_settings))
options_button.grid(row=0, column=0, sticky='NSEW', padx=2, pady=2)
#print(default_settings)

#Create the run button
run_button = tkinter.Button(bottomframe, text="Run!", command=lambda:runprog(default_settings))
run_button.grid(row=0, column=1, sticky='NEWS', padx=2, pady=2)
ent3.bind("<Return>", (lambda event : runprog(default_settings)))

#set up an overall grid in the main window that contains 1 column and one row
mainwindow.grid_columnconfigure(0, weight=1)
mainwindow.grid_rowconfigure(0, weight=1)

mainframe.grid_columnconfigure(0, weight=1)
mainframe.grid_rowconfigure(0, weight=1)

#Within the top half of the window, set up a grid containing 3 columns and 3 rows
#This will control the sizing of the labels, text boxes and buttons
topframe.grid_columnconfigure(0, weight=1)
topframe.grid_columnconfigure(1, weight=10)
topframe.grid_columnconfigure(2, weight=1)

topframe.grid_rowconfigure(0, weight=1)
topframe.grid_rowconfigure(1, weight=1)
topframe.grid_rowconfigure(2, weight=1)

bottomframe.grid_columnconfigure(0, weight=1)
bottomframe.grid_columnconfigure(1, weight=1)
bottomframe.grid_rowconfigure(0, weight=1)


mainwindow.mainloop()
