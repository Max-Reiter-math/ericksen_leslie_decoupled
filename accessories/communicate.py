#!/usr/bin/env python
"""
sends telegram messages to me with updates of the simulation
get telebot by 
pip install pyTelegramBotAPI
"""
# importing all required libraries
import telebot
import json

class msg_bot:
    def __init__(self, config_path):
        self.photo_freq = 10
        with open(config_path) as f:
            json_data = json.load(f)
        self.token = json_data["token"]
        self.chat_id = json_data["chat_id"]
        self.bot = telebot.TeleBot(self.token)

    def send_message(self, message):
        self.bot.send_message(self.chat_id, message)

    def send_photo(self, img_path):
        self.bot.send_photo(self.chat_id, open(img_path, 'rb'))