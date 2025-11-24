# модуль для работы с архивами
import gzip
# модуль для работы с файловой системой
import os
# модуль для графического интерфейса
import tkinter as tk
from tkinter import filedialog, messagebox
# модуль для расширенных коллекций
from collections import defaultdict
# модуль для работы с изображениями
from PIL import Image, ImageTk
# модули для загрузки изображения из URL
from urllib.request import urlopen
from io import BytesIO

class FASTQProcessor:
    """класс для обработки fastq файлов"""
    
    def __init__(self):
        # инициализация счетчиков и хранилищ данных
        self.sequences_count = 0          # количество последовательностей
        self.total_length = 0             # общая длина всех последовательностей
        self.sequence_lengths = []        # список длин каждой последовательности
        self.quality_scores = []          # список оценок качества
        self.base_counts = defaultdict(int) # счетчик оснований
        
    def parse_fastq(self, file_path):
        """
        парсинг fastq файла и сбор статистики
        
        параметры:
            file_path - путь к файлу для анализа
        """
        # сброс предыдущей статистики
        self.sequences_count = 0
        self.total_length = 0
        self.sequence_lengths = []
        self.quality_scores = []
        self.base_counts = defaultdict(int)
        
        # выбор способа открытия файла (обычный или gzip)
        if file_path.endswith('.gz'):
            # открытие gzip файла в текстовом режиме
            opener = gzip.open
        else:
            # открытие обычного текстового файла
            opener = open
        
        # открытие файла с обработкой ошибок кодировки
        with opener(file_path, 'rt', encoding='utf-8', errors='ignore') as file:
            # чтение всех строк файла
            lines = file.readlines()
            
        # обработка файла по блокам по 4 строки (стандартный формат fastq)
        for i in range(0, len(lines), 4):
            # проверка что есть достаточно строк для полного блока
            if i + 3 < len(lines):
                # строка с идентификатором (игнорируем)
                # identifier = lines[i].strip()
                
                # строка с последовательностью (вторая строка в блоке)
                sequence_line = lines[i + 1].strip()
                # строка с разделителем (игнорируем)
                # separator = lines[i + 2].strip()
                # строка с качеством (четвертая строка в блоке)
                quality_line = lines[i + 3].strip()
                
                # увеличение счетчика последовательностей
                self.sequences_count += 1
                # вычисление длины текущей последовательности
                seq_len = len(sequence_line)
                # добавление к общей длине
                self.total_length += seq_len
                # сохранение длины последовательности
                self.sequence_lengths.append(seq_len)
                
                # подсчет оснований в последовательности
                self._count_bases(sequence_line)
                # анализ качества последовательности
                self._analyze_quality(quality_line)
    
    def _count_bases(self, sequence):
        """
        подсчет количества каждого типа оснований в последовательности
        
        параметры:
            sequence - строка с последовательностью нуклеотидов
        """
        # преобразование к верхнему регистру для единообразия
        sequence_upper = sequence.upper()
        
        # подсчет каждого основания в последовательности
        for base in sequence_upper:
            if base in 'ATGC':
                self.base_counts[base] += 1
    
    def _analyze_quality(self, quality_string):
        """
        анализ строки с оценками качества
        
        параметры:
            quality_string - строка с символами качества
        """
        # преобразование символов качества в числовые значения
        # в формате fastq качество кодируется как ord(char) - 33
        quality_scores = [ord(char) - 33 for char in quality_string]
        # добавление оценок качества в общий список
        self.quality_scores.extend(quality_scores)
    
    def get_basic_stats(self):
        """получение базовой статистики по файлу"""
        # проверка наличия данных
        if self.sequences_count == 0:
            return "нет данных для анализа"
            
        # вычисление средней длины последовательности
        avg_length = self.total_length / self.sequences_count
        
        # формирование строки со статистикой
        stats = f"общее количество последовательностей: {self.sequences_count}\n"
        stats += f"средняя длина последовательности: {avg_length:.1f} bp\n"
        stats += f"минимальная длина: {min(self.sequence_lengths)} bp\n"
        stats += f"максимальная длина: {max(self.sequence_lengths)} bp\n"
        stats += f"общая длина всех последовательностей: {self.total_length} bp\n"
        
        return stats
    
    def get_quality_stats(self):
        """получение статистики по качеству последовательностей"""
        if not self.quality_scores:
            return "нет данных о качестве"
            
        # вычисление статистик качества
        avg_quality = sum(self.quality_scores) / len(self.quality_scores)
        min_quality = min(self.quality_scores)
        max_quality = max(self.quality_scores)
        
        stats = f"среднее качество: {avg_quality:.2f}\n"
        stats += f"минимальное качество: {min_quality}\n"
        stats += f"максимальное качество: {max_quality}\n"
        stats += f"общее количество оценок качества: {len(self.quality_scores)}"
        
        return stats
    
    def get_base_composition(self):
        """получение состава оснований"""
        if not self.base_counts:
            return "нет данных о составе оснований"
            
        # вычисление общего количества оснований
        total_bases = sum(self.base_counts.values())
        
        # формирование строки с составом оснований
        composition = "состав оснований:\n"
        for base in 'ATGC':
            count = self.base_counts[base]
            percentage = (count / total_bases) * 100 if total_bases > 0 else 0
            composition += f"  {base}: {count} ({percentage:.1f}%)\n"
        
        return composition

class FASTQAnalyzer:
    """класс графического интерфейса для анализа fastq файлов"""
    
    def __init__(self, root):
        """
        инициализация интерфейса
        
        параметры:
            root - корневое окно приложения
        """
        self.root = root
        # установка заголовка окна
        self.root.title("Анализ FASTQ")
        # установка начального размера окна
        self.root.geometry("700x500")
        # разрешение изменения размера окна
        self.root.resizable(True, True)
        
        # создание обработчика fastq файлов
        self.processor = FASTQProcessor()
        # переменная для хранения пути к текущему файлу
        self.current_file = None
        
        # загрузка и установка фонового изображения
        self.setup_background()
        
        # создание элементов интерфейса
        self.create_interface()
        
        # привязка события изменения размера окна
        self.root.bind("<Configure>", self.on_resize)
    
    def setup_background(self):
        """настройка фонового изображения"""
        try:
            # URL изображения
            image_url = "https://i.pinimg.com/736x/61/aa/70/61aa7092cbad2ac980efb71deead080a.jpg"
            
            # загрузка изображения из URL с помощью urllib
            with urlopen(image_url) as response:
                image_data = response.read()
            
            # открытие изображения с помощью PIL
            self.original_image = Image.open(BytesIO(image_data))
            
            # создание canvas для фона
            self.canvas = tk.Canvas(self.root)
            self.canvas.pack(fill="both", expand=True)
            
            # начальная установка фонового изображения
            self.resize_background()
            
        except Exception as e:
            print(f"Ошибка загрузки фонового изображения: {e}")
            # создание обычного canvas если изображение не загружено
            self.canvas = tk.Canvas(self.root, width=700, height=500, bg="lightgray")
            self.canvas.pack(fill="both", expand=True)
    
    def resize_background(self):
        """изменение размера фонового изображения под размер окна"""
        if hasattr(self, 'original_image'):
            # получение текущего размера окна
            width = self.root.winfo_width()
            height = self.root.winfo_height()
            
            # если окно еще не отобразилось, используем размер по умолчанию
            if width < 10 or height < 10:
                width, height = 700, 500
            
            # изменение размера изображения с сохранением пропорций
            resized_image = self.original_image.resize((width, height), Image.LANCZOS)
            
            # преобразование для tkinter
            self.bg_image = ImageTk.PhotoImage(resized_image)
            
            # обновление фонового изображения на canvas
            self.canvas.delete("background")
            self.canvas.create_image(0, 0, image=self.bg_image, anchor="nw", tags="background")
    
    def on_resize(self, event):
        """обработчик изменения размера окна"""
        # игнорируем события от дочерних виджетов
        if event.widget == self.root:
            self.resize_background()
    
    def create_interface(self):
        """создание пользовательского интерфейса"""
        # создание основного фрейма с отступами поверх canvas
        main_frame = tk.Frame(self.canvas, bg="white", relief="raised", bd=2)
        main_frame.place(relx=0.5, rely=0.5, anchor="center", width=650, height=450)
        
        # заголовок приложения
        title_label = tk.Label(
            main_frame, 
            text="Анализатор файлов FASTQ", 
            font=("DejaVu", 16, "bold"),
            fg="#2c3e50",  # темно-синий цвет в стиле изображения
            bg="white"
        )
        title_label.pack(pady=10)
        
        # фрейм для кнопок управления
        button_frame = tk.Frame(main_frame, bg="white")
        button_frame.pack(pady=10)
        
        # кнопка выбора файла с цветами в стиле изображения
        self.select_btn = tk.Button(
            button_frame, 
            text="выбрать FASTQ файл", 
            command=self.select_file, 
            width=20, 
            height=2,
            bg="#3498db",  # синий цвет в стиле изображения
            fg="white",    # белый текст для контраста
            font=("DejaVu", 10, "bold"),
            relief="raised",
            bd=2,
            activebackground="#2980b9",  # более темный синий при наведении
            activeforeground="white"
        )
        self.select_btn.pack(side=tk.LEFT, padx=5)
        
        # кнопка анализа файла (изначально неактивна) с цветами в стиле изображения
        self.analyze_btn = tk.Button(
            button_frame, 
            text="проанализировать", 
            command=self.analyze_file, 
            state=tk.DISABLED,
            width=15, 
            height=2,
            bg="#2ecc71",  # зеленый цвет для контраста с синим
            fg="white",    # белый текст для контраста
            font=("DejaVu", 10, "bold"),
            relief="raised",
            bd=2,
            activebackground="#27ae60",  # более темный зеленый при наведении
            activeforeground="white"
        )
        self.analyze_btn.pack(side=tk.LEFT, padx=5)
        
        # метка для отображения информации о выбранном файле
        self.file_label = tk.Label(
            main_frame, 
            text="файл не выбран", 
            wraplength=600, 
            justify="left",
            bg="#ecf0f1",  # светло-серый фон
            fg="#2c3e50",  # темно-синий текст
            relief="sunken",
            padx=10,
            pady=5,
            font=("DejaVu", 9)
        )
        self.file_label.pack(pady=5, fill=tk.X, padx=20)
        
        # создание текстовых областей для разных типов статистики
        self.create_text_areas(main_frame)
    
    def create_text_areas(self, parent):
        """создание текстовых областей для отображения статистики"""
        # фрейм для текстовых областей
        text_frame = tk.Frame(parent, bg="white")
        text_frame.pack(pady=10, fill=tk.BOTH, expand=True, padx=20)
        
        # создание текстовых виджетов с метками
        self.create_text_widget(text_frame, "базовая статистика", 0)
        self.create_text_widget(text_frame, "качество последовательностей", 1)
        self.create_text_widget(text_frame, "состав оснований", 2)
    
    def create_text_widget(self, parent, title, row):
        """
        создание текстового виджета с заголовком
        
        параметры:
            parent - родительский фрейм
            title - заголовок текстовой области
            row - номер строки для размещения
        """
        # создание фрейма для текстового виджета
        frame = tk.Frame(parent, bg="white")
        frame.grid(row=row, column=0, sticky="ew", pady=5)
        frame.columnconfigure(0, weight=1)
        
        # метка с заголовком
        label = tk.Label(
            frame, 
            text=title, 
            font=("DejaVu", 11, "bold"), 
            bg="white",
            fg="#2c3e50"  # темно-синий цвет в стиле изображения
        )
        label.grid(row=0, column=0, sticky="w", pady=(0, 5))
        
        # текстовое поле для отображения статистики
        text_widget = tk.Text(
            frame, 
            height=4, 
            width=80, 
            font=("Courier", 9),
            wrap=tk.WORD,
            bg="#f8f9fa",  # очень светлый серый фон
            fg="#2c3e50",  # темно-синий текст
            relief="sunken",
            bd=1,
            selectbackground="#3498db",  # синий цвет выделения
            selectforeground="white"     # белый текст выделения
        )
        text_widget.grid(row=1, column=0, sticky="ew", padx=5)
        
        # добавление скроллбара для текстового поля
        scrollbar = tk.Scrollbar(
            frame, 
            command=text_widget.yview,
            bg="#bdc3c7"  # серый цвет скроллбара
        )
        scrollbar.grid(row=1, column=1, sticky="ns")
        text_widget.config(yscrollcommand=scrollbar.set)
        
        # сохранение ссылки на текстовый виджет
        if row == 0:
            self.basic_text = text_widget
        elif row == 1:
            self.quality_text = text_widget
        elif row == 2:
            self.composition_text = text_widget
    
    def select_file(self):
        """выбор файла через диалоговое окно"""
        # определение фильтров для типов файлов
        file_types = [
            ("fastq files", "*.fastq *.fq"),
            ("gzipped fastq", "*.fastq.gz *.fq.gz"),
            ("all files", "*.*")
        ]
        
        # открытие диалогового окна выбора файла
        filename = filedialog.askopenfilename(
            title="выберите FASTQ файл",
            filetypes=file_types
        )
        
        # обработка выбранного файла
        if filename:
            self.current_file = filename
            # обновление метки с информацией о файле
            file_info = f"выбран файл: {os.path.basename(filename)}"
            self.file_label.config(text=file_info)
            # активация кнопки анализа
            self.analyze_btn.config(state=tk.NORMAL)
    
    def analyze_file(self):
        """анализ выбранного файла"""
        # проверка что файл выбран
        if not self.current_file:
            messagebox.showwarning("внимание", "сначала выберите файл")
            return
            
        try:
            # отображение сообщения о начале обработки
            self.show_loading_message()
            
            # анализ файла с помощью процессора
            self.processor.parse_fastq(self.current_file)
            
            # отображение результатов анализа
            self.display_results()
            
        except Exception as e:
            # обработка ошибок при анализе
            error_message = f"ошибка при обработке файла:\n{str(e)}"
            messagebox.showerror("ошибка", error_message)
            self.clear_results()
    
    def show_loading_message(self):
        """отображение сообщения о загрузке"""
        loading_text = "обработка файла...\nпожалуйста, подождите."
        
        # очистка всех текстовых полей и отображение сообщения о загрузке
        for text_widget in [self.basic_text, self.quality_text, self.composition_text]:
            text_widget.delete(1.0, tk.END)
            text_widget.insert(tk.END, loading_text)
        
        # обновление интерфейса
        self.root.update()
    
    def display_results(self):
        """отображение результатов анализа"""
        # получение статистики из процессора
        basic_stats = self.processor.get_basic_stats()
        quality_stats = self.processor.get_quality_stats()
        composition_stats = self.processor.get_base_composition()
        
        # получение информации о размере файла
        file_size = os.path.getsize(self.current_file)
        file_info = f"\n\nинформация о файле:\n"
        file_info += f"размер файла: {file_size / 1024 / 1024:.2f} MB\n"
        file_info += f"имя файла: {os.path.basename(self.current_file)}"
        
        # отображение базовой статистики
        self.basic_text.delete(1.0, tk.END)
        self.basic_text.insert(tk.END, basic_stats + file_info)
        
        # отображение статистики качества
        self.quality_text.delete(1.0, tk.END)
        self.quality_text.insert(tk.END, quality_stats)
        
        # отображение состава оснований
        self.composition_text.delete(1.0, tk.END)
        self.composition_text.insert(tk.END, composition_stats)
    
    def clear_results(self):
        """очистка результатов анализа"""
        # очистка всех текстовых полей
        self.basic_text.delete(1.0, tk.END)
        self.quality_text.delete(1.0, tk.END)
        self.composition_text.delete(1.0, tk.END)

def main():
    """основная функция запуска приложения"""
    # создание корневого окна
    root = tk.Tk()
    # создание экземпляра приложения
    app = FASTQAnalyzer(root)
    # запуск главного цикла обработки событий
    root.mainloop()

# точка входа в программу
if __name__ == "__main__":

    main()
