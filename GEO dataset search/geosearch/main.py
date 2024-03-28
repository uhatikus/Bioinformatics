import time
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from chromedriver_py import binary_path
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from bs4 import BeautifulSoup

class GEOsearch:
    def __init__(self):
        self.url = "https://artyomovlab.wustl.edu/phantasus/"


    def init_driver(self):
        chrome_options = webdriver.ChromeOptions()
        # chrome_options.add_argument("--headless")
        # chrome_options.add_argument("--window-size=1920,1080")
        # chrome_options.add_argument("--disable-gpu")
        service = Service(executable_path=binary_path)
        driver = webdriver.Chrome(service=service, options=chrome_options)
        return driver

    def parse_website(self):
        driver = self.init_driver()
        driver.get(self.url)
        try:
            time.sleep(3)
            html = driver.page_source
            return html
        except Exception as e:
            print(f"An error occurred: {e}")
        finally:
            driver.quit()
            
    def choose_geo_datasets(self, datasets_to_test):
        good_GEO_datasets = []
        
        for i, GEOdataset in enumerate(datasets_to_test):
            print(f"{int((100*i)/len(datasets_to_test))}% done")
            driver = self.init_driver()
            wait = WebDriverWait(driver, 10)  # Adjust the timeout as needed
            try:
                driver.get(self.url)
                # Wait for the dropdown to be clickable
                # Click the dropdown button
                dropdown_button = wait.until(EC.element_to_be_clickable((By.XPATH, '//button[@title="Choose a file..."]')))
                dropdown_button.click()
                time.sleep(0.5)
                # Select "GEO Datasets" from the dropdown menu
                geo_datasets_option = wait.until(EC.element_to_be_clickable((By.XPATH, '//span[text()="GEO Datasets"]')))
                geo_datasets_option.click()
                time.sleep(0.5)
                # Wait for the input field to appear
                file_geo_input = wait.until(EC.visibility_of_element_located((By.XPATH, '//div[@id="file_geo"]/input[@name="file_geo"]')))
                # Input your GSE or GDS identifier
                file_geo_input.send_keys(GEOdataset)
                # Click the Load button
                load_button = wait.until(EC.element_to_be_clickable((By.XPATH, '//div[@id="file_geo"]/input[@value="Load"]')))
                load_button.click()
                time.sleep(10)
                html = driver.page_source
                # Parse HTML using BeautifulSoup
                soup = BeautifulSoup(html, 'html.parser')
                # Find the div with id "vis"
                vis_div = soup.find('div', id='vis')
                # Find the ul element inside the div with class "nav-tabs"
                ul_element = vis_div.find('ul', class_='nav-tabs')
                # Find the li element with class "phantasus-sortable active"
                li_element = ul_element.find('li', class_='phantasus-sortable active')
                # Get the title attribute of the 'a' element inside the li element
                title = li_element.find('a')['title']
                # Extract the number of rows from the title
                num_rows = int(title.split(' ')[0])
                print(num_rows)
                if num_rows > 10:
                    print(GEOdataset)
                    good_GEO_datasets.append(GEOdataset) 
                # return html
            except Exception as e:
                print(f"An error with {GEOdataset} occurred: {e}")
            finally:
                driver.quit()
        print("100% done")
        return good_GEO_datasets
    
    def save_to_html(self, current_results, filename="output.html"):
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(current_results)


if __name__ == "__main__":
    parser = GEOsearch()
    # datasets_to_test = ["GSE53986", "GSE246672"]
    datasets_to_test = ["GSE248600" ,"GSE248599" ,"GSE241926" ,"GSE241925" ,"GSE241924" ,"GSE227974" ,"GSE227973" ,"GSE227972" ,"GSE227612" ,"GSE218716" ,"GSE218715" ,"GSE260915" ,"GSE260575" ,"GSE193056" ,"GSE188242" ,"GSE260598" ,"GSE260597" ,"GSE260487" ,"GSE246045" ,"GSE246044" ,"GSE241867" ,"GSE241738" ,"GSE240505" ,"GSE240504" ,"GSE222265" ,"GSE222264" ,"GSE207360" ,"GSE197569" ,"GSE164455" ,"GSE164454" ,"GSE164451" ,"GSE164450" ,"GSE260674" ,"GSE259344" ,"GSE226733" ,"GSE226548" ,"GSE211425" ,"GSE168121" ,"GSE260697" ,"GSE249470" ,"GSE259311" ,"GSE256462" ,"GSE256292" ,"GSE252935" ,"GSE252934" ,"GSE252933" ,"GSE249909" ,"GSE249130" ,"GSE248896" ,"GSE247312" ,"GSE240747" ,"GSE239991" ,"GSE238257" ,"GSE235534" ,"GSE233479" ,"GSE230526" ,"GSE230466" ,"GSE230293" ,"GSE230292" ,"GSE230277" ,"GSE230087" ,"GSE228656" ,"GSE228553" ,"GSE228375" ,"GSE227940" ,"GSE227707" ,"GSE227705" ,"GSE227704" ,"GSE226519" ,"GSE221906" ,"GSE217204" ,"GSE213504" ,"GSE200334" ,"GSE198142" ,"GSE197901" ,"GSE197769" ,"GSE255415" ,"GSE255414" ,"GSE255413" ,"GSE255412" ,"GSE255411" ,"GSE255410" ,"GSE255408" ,"GSE255386" ,"GSE255353" ,"GSE255154" ,"GSE255000" ,"GSE253795" ,"GSE249053" ,"GSE243375" ,"GSE218249" ,"GSE206240" ,"GSE256457" ,"GSE245674" ,"GSE241847" ,"GSE240671" ,"GSE236847" ,"GSE229509" ,"GSE226789" ,"GSE254685" ,"GSE254684" ,"GSE246405" ,"GSE227427" ,"GSE201655" ,"GSE256294" ,"GSE249049" ,"GSE248835" ,"GSE247146" ,"GSE245770" ,"GSE245769" ,"GSE245765" ,"GSE244455" ,"GSE244198" ,"GSE244197" ,"GSE243617" ,"GSE237662" ,"GSE237553" ,"GSE237144" ,"GSE236897" ,"GSE236896" ,"GSE234424" ,"GSE227666" ,"GSE226059" ,"GSE220728" ,"GSE255308" ,"GSE224500" ,"GSE220143" ,"GSE167512" ,"GSE256086" ,"GSE256074" ,"GSE256103" ,"GSE256100" ,"GSE256094" ,"GSE251644" ,"GSE241874" ,"GSE233609" ,"GSE225900" ,"GSE254701" ,"GSE254479" ,"GSE249914" ,"GSE249587" ,"GSE244311" ,"GSE235329" ,"GSE167321" ,"GSE256235" ,"GSE248874" ,"GSE244042" ,"GSE240339" ,"GSE237403" ,"GSE228782" ,"GSE228781" ,"GSE256127" ,"GSE256012" ,"GSE256010" ,"GSE256003" ,"GSE255958" ,"GSE255896" ,"GSE255853" ,"GSE255852" ,"GSE255850" ,"GSE255629" ,"GSE255615" ,"GSE255598" ,"GSE244471" ,"GSE244470" ,"GSE232007" ,"GSE225195" ,"GSE211565" ,"GSE256061" ,"GSE256060" ,"GSE256049" ,"GSE242778" ,"GSE238198" ,"GSE238197" ,"GSE227976" ,"GSE227756" ,"GSE224781" ,"GSE224465" ,"GSE197798" ,"GSE255729" ,"GSE255707" ,"GSE255706" ,"GSE244327" ,"GSE244326" ,"GSE244287" ,"GSE242717" ,"GSE255752" ,"GSE255143" ,"GSE231560" ,"GSE228430" ,"GSE255970" ,"GSE255894" ,"GSE255187" ,"GSE252645" ,"GSE250348" ,"GSE246091" ,"GSE241422" ,"GSE240055" ,"GSE222315" ,"GSE201536" ,"GSE255898" ,"GSE252155" ,"GSE251743" ,"GSE249098" ,"GSE247951" ,"GSE247949" ,"GSE247948" ,"GSE247947" ,"GSE241895" ,"GSE241883" ,"GSE235648" ,"GSE229877" ,"GSE221062" ,"GSE221056" ,"GSE216810" ,"GSE206336" ,"GSE255448" ,"GSE255286" ,"GSE255129" ,"GSE254035" ,"GSE252806" ,"GSE252406" ,"GSE249569" ,"GSE248011" ,"GSE247460" ,"GSE241871" ,"GSE239474" ,"GSE238207" ,"GSE238108" ,"GSE236505" ,"GSE236504" ,"GSE236503" ,"GSE234866" ,"GSE225210" ,"GSE184072" ,"GSE153033" ,"GSE255581" ,"GSE253400" ,"GSE244353" ,"GSE241224" ,"GSE216117" ,"GSE213762" ,"GSE199503" ,"GSE255589" ,"GSE255313" ,"GSE255080" ,"GSE249536" ,"GSE249316" ,"GSE236026" ,"GSE234475" ,"GSE230362" ,"GSE229666" ,"GSE229664" ,"GSE229663" ,"GSE225352" ,"GSE202397" ,"GSE181905" ,"GSE255300" ,"GSE244626" ,"GSE196689" ,"GSE255107" ,"GSE254649" ,"GSE249825" ,"GSE235908" ,"GSE235555" ,"GSE240305" ,"GSE235549" ,"GSE235031" ,"GSE232609" ,"GSE225022" ,"GSE225021" ,"GSE225017" ,"GSE225013" ,"GSE225006" ,"GSE255013" ,"GSE242154" ,"GSE241922" ,"GSE254937" ,"GSE254718" ,"GSE245011" ,"GSE240722" ,"GSE239295" ,"GSE239294" ,"GSE238116" ,"GSE237033" ,"GSE235621" ,"GSE235374" ,"GSE228680" ,"GSE227907" ,"GSE220879" ,"GSE190738" ,"GSE255165" ,"GSE255163" ,"GSE255058" ,"GSE254833" ,"GSE254832" ,"GSE254829" ,"GSE252051" ,"GSE250050" ,"GSE244983" ,"GSE244982" ,"GSE243568" ,"GSE243565" ,"GSE243564" ,"GSE243562" ,"GSE232572" ,"GSE232571" ,"GSE254906" ,"GSE254905" ,"GSE254904" ,"GSE254902" ,"GSE254900" ,"GSE254899" ,"GSE254898" ,"GSE254896" ,"GSE254851" ,"GSE254526"]
    results = parser.choose_geo_datasets(datasets_to_test)
    print(results)
    # if results:
    #     parser.save_to_html(results, filename="phantasus.html")
    #     print("Website contents saved to phantasus.html")
    # else:
    #     print("No results to save.")
    with open("results.txt", 'w') as f:
        f.writelines(f"{item}\n" for item in results)