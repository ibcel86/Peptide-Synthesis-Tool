import customtkinter

class TabView(customtkinter.CTkTabview):
    def __init__(self, master):
        super().__init__(master)

        # Tabs
        self.add("Synthesis Planner")
        self.add("Modify Synthesis")

        self.title_synthesisplanner = customtkinter.CTkLabel(self.tab("Synthesis Planner"), text="Synthesis Planner")
        self.title_synthesisplanner.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        self.title_modifysynthesis = customtkinter.CTkLabel(self.tab("Modify Synthesis"),text="Modify Synthesis")
        self.title_modifysynthesis.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="e")


class App(customtkinter.CTk):
    '''Activates the GUI etc'''
    def __init__(self):
        super().__init__()

        self.title("Protein Sequence Calculator")
        customtkinter.set_appearance_mode("dark")

        self.geometry("1280x720")
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)

        self.tabview = TabView(self)
        self.tabview.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

