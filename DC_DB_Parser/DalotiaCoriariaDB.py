import tkinter as tk
import tkinter.ttk as ttk
import tkinter.filedialog as filedialog
import tkinter.messagebox as messagebox
from dalotia_db.dalotia_db import main
from dalotia_db.dalotia_db import dalotia_db as ddbc
import subprocess
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq

class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.master.title("Dalotia Coriaria Probe Generator")
        self.pack()
        self.create_widgets()

    def create_widgets(self):
        self.fasta_label = tk.Label(self, text="FASTA File:")
        self.fasta_label.grid(row=0, column=0)
        self.fasta_entry = tk.Entry(self)
        self.fasta_entry.grid(row=0, column=1)

        self.gff_label = tk.Label(self, text="GFF File:")
        self.gff_label.grid(row=1, column=0)
        self.gff_entry = tk.Entry(self)
        self.gff_entry.grid(row=1, column=1)

        self.genelist_label = tk.Label(self, text="Gene List ( as evm.model.hic_scaffold_x_pilon.x):")
        self.genelist_label.grid(row=2, column=0)
        self.genelist_entry = tk.Entry(self)
        self.genelist_entry.grid(row=2, column=1)

        self.biotype_label = tk.Label(self, text="Biotype:")
        self.biotype_label.grid(row=3, column=0)
        self.biotype_var = tk.StringVar(value='mRNA')
        self.biotype_dropdown = ttk.Combobox(self, textvariable=self.biotype_var, values=['mRNA', 'exons', 'introns'])
        self.biotype_dropdown.grid(row=3, column=1)
        self.biotype_dropdown.bind("<<ComboboxSelected>>", self.enable_bed_entry)

        self.bed_label = tk.Label(self, text="BED File (optional):")
        self.bed_label.grid(row=4, column=0)
        self.bed_entry = tk.Entry(self, state=tk.DISABLED)
        self.bed_entry.grid(row=4, column=1)

        self.probes_label = tk.Label(self, text="Probes CSV(after probe generation):")
        self.probes_label.grid(row=5, column=0)
        self.probes_entry = tk.Entry(self)
        self.probes_entry.grid(row=5, column=1)

        self.select_fasta_file_button = tk.Button(self, text="Select FASTA File", command=self.select_fasta_file)
        self.select_fasta_file_button.grid(row=0, column=2)

        self.select_gff_file_button = tk.Button(self, text="Select GFF File", command=self.select_gff_file)
        self.select_gff_file_button.grid(row=1, column=2)

        self.select_genelist_file_button = tk.Button(self, text="Select Gene List File", command=self.select_genelist_file)
        self.select_genelist_file_button.grid(row=2, column=2)

        self.select_genelist_file_button = tk.Button(self, text="Select Probes List File", command=self.select_probes_file)
        self.select_genelist_file_button.grid(row=5, column=2)

        self.generate_button = tk.Button(self, text="Generate FASTA", command=self.generate_fasta)
        self.generate_button.grid(row=7, column=0)

        self.generate_button = tk.Button(self, text="BLAST", command=self.blast_probes)
        self.generate_button.grid(row=7, column=1)

        self.quit_button = tk.Button(self, text="Quit", fg="red", command=self.master.destroy)
        self.quit_button.grid(row=7, column=2)
    
    def enable_bed_entry(self, *args):
        if self.biotype_var.get() == "mRNA":
            self.bed_entry.config(state=tk.DISABLED)
            self.bed_label.grid_forget()
            self.bed_entry.grid_forget()
        else:
            self.bed_entry.config(state=tk.NORMAL)
            self.bed_label.grid(row=4, column=0)
            self.bed_entry.grid(row=4, column=1)
            
            # ask for BED file if not already selected
            if not self.bed_entry.get():
                messagebox.showinfo("", "Select BED File Location")
                filename = filedialog.askopenfilename()
                self.bed_entry.delete(0, tk.END)
                self.bed_entry.insert(0, filename)

    def select_fasta_file(self):
        filename = filedialog.askopenfilename()
        self.fasta_entry.delete(0, tk.END)
        self.fasta_entry.insert(0, filename)

    def select_gff_file(self):
        filename = filedialog.askopenfilename()
        self.gff_entry.delete(0, tk.END)
        self.gff_entry.insert(0, filename)

    def select_genelist_file(self):
        filename = filedialog.askopenfilename()
        self.genelist_entry.delete(0, tk.END)
        self.genelist_entry.insert(0, filename)

    def select_probes_file(self):
        filename = filedialog.askopenfilename()
        self.probes_entry.delete(0, tk.END)
        self.probes_entry.insert(0, filename)

    def generate_fasta(self):
        fasta_file = self.fasta_entry.get()
        gff_file = self.gff_entry.get()
        genelist_file = self.genelist_entry.get()
        biotype = self.biotype_var.get()
        bed_file = self.bed_entry.get() if biotype.lower() in ['exons', 'introns'] else None
        filename = filedialog.asksaveasfilename(defaultextension='.fasta')

        if filename:
            main(fasta_file, gff_file, biotype, genelist_file, bed_file, filename)
            # show a dialog box with a message when finished
            messagebox.showinfo("Finished", "FASTA file generated successfully!")
    
    def blast_probes(self):
        probes_file = self.probes_entry.get()
        seq_only = pd.read_csv(probes_file)[['<GENE>', '<PROBE>']]
        # add a unique suffix 
        for i in range(len(seq_only)):
            seq_only.at[i, '<GENE>'] = str(seq_only.at[i, '<GENE>']) + '_' + str(i)
        # get reverse complement of probes to blast against cDNA
        rev = []
        for probe in seq_only['<PROBE>']:
            rev.append(str(Seq(probe).reverse_complement()))
        seq_only['<PROBE>'] = rev
        output_dir = str(Path(probes_file).parent / "blastfile.fa")
        ddbc.create_fasta(output_dir, seq_only['<PROBE>'].values.tolist(), 
                          seq_only['<GENE>'].values.tolist(), 
                          seq_only['<GENE>'].values.tolist(), biotype='mRNA')
        output_dir2 = str(Path(probes_file).parent / "results.txt")
        bashCommand = f"blastn -db /home/klcolon/dalotia_blastdb/dalotia_blastdb -query {output_dir} -outfmt 6 -out {output_dir2} -evalue 0.001"
        subprocess.run(bashCommand, shell=True)
        df = pd.read_csv(output_dir2, sep = '\t', header = None)
        offtarget_idx = []
        #look for off target > 16 bp
        for qseqid, sseqid, overlap in df[[0,1, 3]].values:
            if sseqid in qseqid:
                continue
            else:
                if overlap < 17:
                    continue
                else:
                    offtarget_idx.append(int(qseqid.split('_')[-1]))
        final_probes = pd.read_csv(probes_file).drop(offtarget_idx)
        output_dir3 = str(Path(probes_file).parent / "filtered_final_probes.csv")
        final_probes.to_csv(output_dir3 , index=False)
        messagebox.showinfo("Probes are ready for use!", " Look for filtered_final_probes.csv")
root = tk.Tk()
app = Application(master=root)
app.mainloop()

