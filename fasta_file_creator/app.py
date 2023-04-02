import tkinter as tk
import tkinter.ttk as ttk
import tkinter.filedialog as filedialog
import tkinter.messagebox as messagebox
from dalotia_db.dalotia_db import main

class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.master.title("Fasta File Generator")
        self.pack()
        self.create_widgets()

    def create_widgets(self):
        self.fasta_label = tk.Label(self, text="Fasta File:")
        self.fasta_label.grid(row=0, column=0)
        self.fasta_entry = tk.Entry(self)
        self.fasta_entry.grid(row=0, column=1)

        self.gff_label = tk.Label(self, text="GFF File:")
        self.gff_label.grid(row=1, column=0)
        self.gff_entry = tk.Entry(self)
        self.gff_entry.grid(row=1, column=1)

        self.biotype_label = tk.Label(self, text="Biotype:")
        self.biotype_label.grid(row=2, column=0)
        self.biotype_var = tk.StringVar(value='mRNA')
        self.biotype_dropdown = ttk.Combobox(self, textvariable=self.biotype_var, values=['mRNA', 'exons', 'introns'])
        self.biotype_dropdown.grid(row=2, column=1)
        self.biotype_dropdown.bind("<<ComboboxSelected>>", self.enable_bed_entry)

        if self.biotype_var.get() == "mRNA":
            self.bed_label = tk.Label(self, text="BED File (optional):")
            self.bed_label.grid(row=3, column=0)
            self.bed_entry = tk.Entry(self, state=tk.DISABLED)
            self.bed_entry.grid(row=3, column=1)

        self.select_fasta_file_button = tk.Button(self, text="Select Fasta File", command=self.select_fasta_file)
        self.select_fasta_file_button.grid(row=0, column=2)

        self.select_gff_file_button = tk.Button(self, text="Select GFF File", command=self.select_gff_file)
        self.select_gff_file_button.grid(row=1, column=2)

        self.generate_button = tk.Button(self, text="Generate Fasta", command=self.generate_fasta)
        self.generate_button.grid(row=4, column=1)

        self.quit_button = tk.Button(self, text="Quit", fg="red",
                            command=self.master.destroy)
        self.quit_button.grid(row=4, column=2)

    def enable_bed_entry(self, *args):
        if self.biotype_var.get() in ["introns", "exons"]:
            self.bed_entry.config(state=tk.NORMAL)
        else:
            self.bed_entry.delete(0, tk.END)
            self.bed_entry.config(state=tk.DISABLED)

    def select_fasta_file(self):
        filename = filedialog.askopenfilename()
        self.fasta_entry.delete(0, tk.END)
        self.fasta_entry.insert(0, filename)

    def select_gff_file(self):
        filename = filedialog.askopenfilename()
        self.gff_entry.delete(0, tk.END)
        self.gff_entry.insert(0, filename)

    def generate_fasta(self):
        fasta_file = self.fasta_entry.get()
        gff_file = self.gff_entry.get()
        biotype = self.biotype_var.get()
        bed_file = self.bed_entry.get() if biotype.lower() in ['exons', 'introns'] else None

        main(fasta_file, gff_file, biotype, 'output.fasta', bed_file)
         # show a dialog box with a message when finished
        messagebox.showinfo("Finished", "Fasta file generated successfully!")

root = tk.Tk()
app = Application(master=root)
app.mainloop()

