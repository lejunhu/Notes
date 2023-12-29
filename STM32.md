This is a note of STM32 study
%coded in GBK
# Ŀ¼
- [Ŀ¼](#Ŀ¼)
- [��һ�£�STM32���](#��һ��stm32���)
  - [����](#����)
  - [STM32F103C8T6](#stm32f103c8t6)
    - [���](#���)
    - [Peripheral](#peripheral)
    - [ϵͳ�ṹ](#ϵͳ�ṹ)
    - [���Ŷ���](#���Ŷ���)

# ��һ�£�STM32���

## ����

STM32��ST��˾���� ARM Cortex-M �ں˿�����32λ΢��������

����Cortexϵ���ںˣ���

## STM32F103C8T6
### ���
- ϵ�У�Mainstream STM32F1
- �ںˣ�ARM Cortex M3
- ��Ƶ��72MHz
- RAM: 20K(SRAM)
- ROM: 64K(Flash)
- ���磺2.0~3.6V(standard 3.3V)
- ��װ��LQFP48
### Peripheral
| Abbreviation |                          Name                          |      Note       |
| :----------: | :----------------------------------------------------: | :-------------: |
|     NVIC     | Nested Vectored Interrupt ControllerǶ�������жϿ����� |    in kernel    |
|   Systick    |                     ϵͳ�δ�ʱ��                     |    in kernel    |
|     RCC      |                     ��λ��ʱ�ӿ���                     |
|     GPIO     |                        ͨ��IO��                        |
|     AFIO     |                        ����IO��                        |
|     EXTI     |                        �ⲿ�ж�                        |
|     TIM      |                         ��ʱ��                         |
|     ADC      |                       ģ��ת����                       |
|     DMA      |                      ֱ���ڴ����                      |
|    USART     |                   ͬ��/�첽����ͨ��                    |
|     I2C      |                        I2C Com                         |
|     SPI      |                        SPIͨ��                         |
|     CAN      |                        CANͨ��                         |
|     USB      |                        USBͨ��                         |
|     RTC      |                        ʵʱʱ��                        |
|     CRC      |                        CRCУ��                         |
|     PWR      |                        ��Դ����                        |
|     BKP      |                       ���ݼĴ���                       |
|     IWDG     |                       �������Ź�                       |
|     WWDG     |                       ���ڿ��Ź�                       |
|     DAC      |                       ��ģת����                       | not in f103c8t6 |
|     SDIO     |                        SD���ӿ�                        | not in f103c8t6 |
|     FSMC     |                   �ɱ侲̬�洢������                   | not in f103c8t6 |
|   USB OTG    |                      USB�����ӿ�                       | not in f103c8t6 |
### ϵͳ�ṹ
![ϵͳ�ṹͼ](STM32F103_Struct.png)
- ICode:ָ�����ߣ���Ҫ����FLASH�����س���ָ��
- DCode:�������ߣ���Ҫ����FLASH����������
- System:ϵͳ���ߣ�����SRAM,FSMC����������
- AHBϵͳ���ߣ����ڹ�����Ҫ����
### ���Ŷ���
![���Ŷ���](pins.png)
### ��Сϵͳ��·
��

# �ڶ���:GPIO
��

# �����£�EXTI�ⲿ�ж�
## ���
EXTI��APB2�����ϵ����裬���Լ��ָ��GPIO�ڵĵ�ƽ�źţ�����ָ����GPIO�ڲ�����ƽ�仯ʱ��EXTI��������NVIC�����ж����룬����NVIC�þ�������ж�������ִ�ж�Ӧ��EXTI�жϳ���
