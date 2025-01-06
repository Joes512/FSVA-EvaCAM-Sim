import random

def extract_content(base, length, output_file):
    try:
        # 讀取 reference.txt
        with open('reference.txt', 'r', encoding='utf-8') as ref_file:
            content = ref_file.read()

        # 確保 base 和 length 不超過內容長度
        if base < 0 or base >= len(content):
            raise ValueError("Base index is out of bounds.")

        if length < 0:
            raise ValueError("Length must be a non-negative integer.")

        extracted = content[base:base + length]

        # 將擷取的內容寫入指定檔案
        with open(output_file, 'w', encoding='utf-8') as read_file:
            read_file.write(extracted)

        print(f"Content extracted successfully and saved to {output_file}.")

    except FileNotFoundError:
        print("Error: 'reference.txt' not found in the current directory.")
    except ValueError as ve:
        print(f"Error: {ve}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    try:
        # 隨機產生 10 個 base 和對應的長度
        random_bases = [random.randint(0, 1000) for _ in range(10)]
        lengths = [10000] * 10 + [200] * 10

        for i, base in enumerate(random_bases * 2):
            length = lengths[i]
            output_file = f"{base}_{length}.txt"
            extract_content(base, length, output_file)

    except ValueError:
        print("An error occurred during the extraction process.")
