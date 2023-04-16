
build/test/TEST_gather_scatter:     file format elf64-x86-64


Disassembly of section .init:

0000000000009000 <_init>:
    9000:	f3 0f 1e fa          	endbr64 
    9004:	48 83 ec 08          	sub    $0x8,%rsp
    9008:	48 8b 05 d9 3f 01 00 	mov    0x13fd9(%rip),%rax        # 1cfe8 <__gmon_start__@Base>
    900f:	48 85 c0             	test   %rax,%rax
    9012:	74 02                	je     9016 <_init+0x16>
    9014:	ff d0                	call   *%rax
    9016:	48 83 c4 08          	add    $0x8,%rsp
    901a:	c3                   	ret    

Disassembly of section .plt:

0000000000009020 <.plt>:
    9020:	ff 35 ca 38 01 00    	push   0x138ca(%rip)        # 1c8f0 <_GLOBAL_OFFSET_TABLE_+0x8>
    9026:	f2 ff 25 cb 38 01 00 	bnd jmp *0x138cb(%rip)        # 1c8f8 <_GLOBAL_OFFSET_TABLE_+0x10>
    902d:	0f 1f 00             	nopl   (%rax)
    9030:	f3 0f 1e fa          	endbr64 
    9034:	68 00 00 00 00       	push   $0x0
    9039:	f2 e9 e1 ff ff ff    	bnd jmp 9020 <_init+0x20>
    903f:	90                   	nop
    9040:	f3 0f 1e fa          	endbr64 
    9044:	68 01 00 00 00       	push   $0x1
    9049:	f2 e9 d1 ff ff ff    	bnd jmp 9020 <_init+0x20>
    904f:	90                   	nop
    9050:	f3 0f 1e fa          	endbr64 
    9054:	68 02 00 00 00       	push   $0x2
    9059:	f2 e9 c1 ff ff ff    	bnd jmp 9020 <_init+0x20>
    905f:	90                   	nop
    9060:	f3 0f 1e fa          	endbr64 
    9064:	68 03 00 00 00       	push   $0x3
    9069:	f2 e9 b1 ff ff ff    	bnd jmp 9020 <_init+0x20>
    906f:	90                   	nop
    9070:	f3 0f 1e fa          	endbr64 
    9074:	68 04 00 00 00       	push   $0x4
    9079:	f2 e9 a1 ff ff ff    	bnd jmp 9020 <_init+0x20>
    907f:	90                   	nop
    9080:	f3 0f 1e fa          	endbr64 
    9084:	68 05 00 00 00       	push   $0x5
    9089:	f2 e9 91 ff ff ff    	bnd jmp 9020 <_init+0x20>
    908f:	90                   	nop
    9090:	f3 0f 1e fa          	endbr64 
    9094:	68 06 00 00 00       	push   $0x6
    9099:	f2 e9 81 ff ff ff    	bnd jmp 9020 <_init+0x20>
    909f:	90                   	nop
    90a0:	f3 0f 1e fa          	endbr64 
    90a4:	68 07 00 00 00       	push   $0x7
    90a9:	f2 e9 71 ff ff ff    	bnd jmp 9020 <_init+0x20>
    90af:	90                   	nop
    90b0:	f3 0f 1e fa          	endbr64 
    90b4:	68 08 00 00 00       	push   $0x8
    90b9:	f2 e9 61 ff ff ff    	bnd jmp 9020 <_init+0x20>
    90bf:	90                   	nop
    90c0:	f3 0f 1e fa          	endbr64 
    90c4:	68 09 00 00 00       	push   $0x9
    90c9:	f2 e9 51 ff ff ff    	bnd jmp 9020 <_init+0x20>
    90cf:	90                   	nop
    90d0:	f3 0f 1e fa          	endbr64 
    90d4:	68 0a 00 00 00       	push   $0xa
    90d9:	f2 e9 41 ff ff ff    	bnd jmp 9020 <_init+0x20>
    90df:	90                   	nop
    90e0:	f3 0f 1e fa          	endbr64 
    90e4:	68 0b 00 00 00       	push   $0xb
    90e9:	f2 e9 31 ff ff ff    	bnd jmp 9020 <_init+0x20>
    90ef:	90                   	nop
    90f0:	f3 0f 1e fa          	endbr64 
    90f4:	68 0c 00 00 00       	push   $0xc
    90f9:	f2 e9 21 ff ff ff    	bnd jmp 9020 <_init+0x20>
    90ff:	90                   	nop
    9100:	f3 0f 1e fa          	endbr64 
    9104:	68 0d 00 00 00       	push   $0xd
    9109:	f2 e9 11 ff ff ff    	bnd jmp 9020 <_init+0x20>
    910f:	90                   	nop
    9110:	f3 0f 1e fa          	endbr64 
    9114:	68 0e 00 00 00       	push   $0xe
    9119:	f2 e9 01 ff ff ff    	bnd jmp 9020 <_init+0x20>
    911f:	90                   	nop
    9120:	f3 0f 1e fa          	endbr64 
    9124:	68 0f 00 00 00       	push   $0xf
    9129:	f2 e9 f1 fe ff ff    	bnd jmp 9020 <_init+0x20>
    912f:	90                   	nop
    9130:	f3 0f 1e fa          	endbr64 
    9134:	68 10 00 00 00       	push   $0x10
    9139:	f2 e9 e1 fe ff ff    	bnd jmp 9020 <_init+0x20>
    913f:	90                   	nop
    9140:	f3 0f 1e fa          	endbr64 
    9144:	68 11 00 00 00       	push   $0x11
    9149:	f2 e9 d1 fe ff ff    	bnd jmp 9020 <_init+0x20>
    914f:	90                   	nop
    9150:	f3 0f 1e fa          	endbr64 
    9154:	68 12 00 00 00       	push   $0x12
    9159:	f2 e9 c1 fe ff ff    	bnd jmp 9020 <_init+0x20>
    915f:	90                   	nop
    9160:	f3 0f 1e fa          	endbr64 
    9164:	68 13 00 00 00       	push   $0x13
    9169:	f2 e9 b1 fe ff ff    	bnd jmp 9020 <_init+0x20>
    916f:	90                   	nop
    9170:	f3 0f 1e fa          	endbr64 
    9174:	68 14 00 00 00       	push   $0x14
    9179:	f2 e9 a1 fe ff ff    	bnd jmp 9020 <_init+0x20>
    917f:	90                   	nop
    9180:	f3 0f 1e fa          	endbr64 
    9184:	68 15 00 00 00       	push   $0x15
    9189:	f2 e9 91 fe ff ff    	bnd jmp 9020 <_init+0x20>
    918f:	90                   	nop
    9190:	f3 0f 1e fa          	endbr64 
    9194:	68 16 00 00 00       	push   $0x16
    9199:	f2 e9 81 fe ff ff    	bnd jmp 9020 <_init+0x20>
    919f:	90                   	nop
    91a0:	f3 0f 1e fa          	endbr64 
    91a4:	68 17 00 00 00       	push   $0x17
    91a9:	f2 e9 71 fe ff ff    	bnd jmp 9020 <_init+0x20>
    91af:	90                   	nop
    91b0:	f3 0f 1e fa          	endbr64 
    91b4:	68 18 00 00 00       	push   $0x18
    91b9:	f2 e9 61 fe ff ff    	bnd jmp 9020 <_init+0x20>
    91bf:	90                   	nop
    91c0:	f3 0f 1e fa          	endbr64 
    91c4:	68 19 00 00 00       	push   $0x19
    91c9:	f2 e9 51 fe ff ff    	bnd jmp 9020 <_init+0x20>
    91cf:	90                   	nop
    91d0:	f3 0f 1e fa          	endbr64 
    91d4:	68 1a 00 00 00       	push   $0x1a
    91d9:	f2 e9 41 fe ff ff    	bnd jmp 9020 <_init+0x20>
    91df:	90                   	nop
    91e0:	f3 0f 1e fa          	endbr64 
    91e4:	68 1b 00 00 00       	push   $0x1b
    91e9:	f2 e9 31 fe ff ff    	bnd jmp 9020 <_init+0x20>
    91ef:	90                   	nop
    91f0:	f3 0f 1e fa          	endbr64 
    91f4:	68 1c 00 00 00       	push   $0x1c
    91f9:	f2 e9 21 fe ff ff    	bnd jmp 9020 <_init+0x20>
    91ff:	90                   	nop
    9200:	f3 0f 1e fa          	endbr64 
    9204:	68 1d 00 00 00       	push   $0x1d
    9209:	f2 e9 11 fe ff ff    	bnd jmp 9020 <_init+0x20>
    920f:	90                   	nop
    9210:	f3 0f 1e fa          	endbr64 
    9214:	68 1e 00 00 00       	push   $0x1e
    9219:	f2 e9 01 fe ff ff    	bnd jmp 9020 <_init+0x20>
    921f:	90                   	nop
    9220:	f3 0f 1e fa          	endbr64 
    9224:	68 1f 00 00 00       	push   $0x1f
    9229:	f2 e9 f1 fd ff ff    	bnd jmp 9020 <_init+0x20>
    922f:	90                   	nop
    9230:	f3 0f 1e fa          	endbr64 
    9234:	68 20 00 00 00       	push   $0x20
    9239:	f2 e9 e1 fd ff ff    	bnd jmp 9020 <_init+0x20>
    923f:	90                   	nop
    9240:	f3 0f 1e fa          	endbr64 
    9244:	68 21 00 00 00       	push   $0x21
    9249:	f2 e9 d1 fd ff ff    	bnd jmp 9020 <_init+0x20>
    924f:	90                   	nop
    9250:	f3 0f 1e fa          	endbr64 
    9254:	68 22 00 00 00       	push   $0x22
    9259:	f2 e9 c1 fd ff ff    	bnd jmp 9020 <_init+0x20>
    925f:	90                   	nop
    9260:	f3 0f 1e fa          	endbr64 
    9264:	68 23 00 00 00       	push   $0x23
    9269:	f2 e9 b1 fd ff ff    	bnd jmp 9020 <_init+0x20>
    926f:	90                   	nop
    9270:	f3 0f 1e fa          	endbr64 
    9274:	68 24 00 00 00       	push   $0x24
    9279:	f2 e9 a1 fd ff ff    	bnd jmp 9020 <_init+0x20>
    927f:	90                   	nop
    9280:	f3 0f 1e fa          	endbr64 
    9284:	68 25 00 00 00       	push   $0x25
    9289:	f2 e9 91 fd ff ff    	bnd jmp 9020 <_init+0x20>
    928f:	90                   	nop
    9290:	f3 0f 1e fa          	endbr64 
    9294:	68 26 00 00 00       	push   $0x26
    9299:	f2 e9 81 fd ff ff    	bnd jmp 9020 <_init+0x20>
    929f:	90                   	nop
    92a0:	f3 0f 1e fa          	endbr64 
    92a4:	68 27 00 00 00       	push   $0x27
    92a9:	f2 e9 71 fd ff ff    	bnd jmp 9020 <_init+0x20>
    92af:	90                   	nop
    92b0:	f3 0f 1e fa          	endbr64 
    92b4:	68 28 00 00 00       	push   $0x28
    92b9:	f2 e9 61 fd ff ff    	bnd jmp 9020 <_init+0x20>
    92bf:	90                   	nop
    92c0:	f3 0f 1e fa          	endbr64 
    92c4:	68 29 00 00 00       	push   $0x29
    92c9:	f2 e9 51 fd ff ff    	bnd jmp 9020 <_init+0x20>
    92cf:	90                   	nop
    92d0:	f3 0f 1e fa          	endbr64 
    92d4:	68 2a 00 00 00       	push   $0x2a
    92d9:	f2 e9 41 fd ff ff    	bnd jmp 9020 <_init+0x20>
    92df:	90                   	nop
    92e0:	f3 0f 1e fa          	endbr64 
    92e4:	68 2b 00 00 00       	push   $0x2b
    92e9:	f2 e9 31 fd ff ff    	bnd jmp 9020 <_init+0x20>
    92ef:	90                   	nop
    92f0:	f3 0f 1e fa          	endbr64 
    92f4:	68 2c 00 00 00       	push   $0x2c
    92f9:	f2 e9 21 fd ff ff    	bnd jmp 9020 <_init+0x20>
    92ff:	90                   	nop
    9300:	f3 0f 1e fa          	endbr64 
    9304:	68 2d 00 00 00       	push   $0x2d
    9309:	f2 e9 11 fd ff ff    	bnd jmp 9020 <_init+0x20>
    930f:	90                   	nop
    9310:	f3 0f 1e fa          	endbr64 
    9314:	68 2e 00 00 00       	push   $0x2e
    9319:	f2 e9 01 fd ff ff    	bnd jmp 9020 <_init+0x20>
    931f:	90                   	nop
    9320:	f3 0f 1e fa          	endbr64 
    9324:	68 2f 00 00 00       	push   $0x2f
    9329:	f2 e9 f1 fc ff ff    	bnd jmp 9020 <_init+0x20>
    932f:	90                   	nop
    9330:	f3 0f 1e fa          	endbr64 
    9334:	68 30 00 00 00       	push   $0x30
    9339:	f2 e9 e1 fc ff ff    	bnd jmp 9020 <_init+0x20>
    933f:	90                   	nop
    9340:	f3 0f 1e fa          	endbr64 
    9344:	68 31 00 00 00       	push   $0x31
    9349:	f2 e9 d1 fc ff ff    	bnd jmp 9020 <_init+0x20>
    934f:	90                   	nop
    9350:	f3 0f 1e fa          	endbr64 
    9354:	68 32 00 00 00       	push   $0x32
    9359:	f2 e9 c1 fc ff ff    	bnd jmp 9020 <_init+0x20>
    935f:	90                   	nop
    9360:	f3 0f 1e fa          	endbr64 
    9364:	68 33 00 00 00       	push   $0x33
    9369:	f2 e9 b1 fc ff ff    	bnd jmp 9020 <_init+0x20>
    936f:	90                   	nop
    9370:	f3 0f 1e fa          	endbr64 
    9374:	68 34 00 00 00       	push   $0x34
    9379:	f2 e9 a1 fc ff ff    	bnd jmp 9020 <_init+0x20>
    937f:	90                   	nop
    9380:	f3 0f 1e fa          	endbr64 
    9384:	68 35 00 00 00       	push   $0x35
    9389:	f2 e9 91 fc ff ff    	bnd jmp 9020 <_init+0x20>
    938f:	90                   	nop
    9390:	f3 0f 1e fa          	endbr64 
    9394:	68 36 00 00 00       	push   $0x36
    9399:	f2 e9 81 fc ff ff    	bnd jmp 9020 <_init+0x20>
    939f:	90                   	nop
    93a0:	f3 0f 1e fa          	endbr64 
    93a4:	68 37 00 00 00       	push   $0x37
    93a9:	f2 e9 71 fc ff ff    	bnd jmp 9020 <_init+0x20>
    93af:	90                   	nop
    93b0:	f3 0f 1e fa          	endbr64 
    93b4:	68 38 00 00 00       	push   $0x38
    93b9:	f2 e9 61 fc ff ff    	bnd jmp 9020 <_init+0x20>
    93bf:	90                   	nop
    93c0:	f3 0f 1e fa          	endbr64 
    93c4:	68 39 00 00 00       	push   $0x39
    93c9:	f2 e9 51 fc ff ff    	bnd jmp 9020 <_init+0x20>
    93cf:	90                   	nop
    93d0:	f3 0f 1e fa          	endbr64 
    93d4:	68 3a 00 00 00       	push   $0x3a
    93d9:	f2 e9 41 fc ff ff    	bnd jmp 9020 <_init+0x20>
    93df:	90                   	nop
    93e0:	f3 0f 1e fa          	endbr64 
    93e4:	68 3b 00 00 00       	push   $0x3b
    93e9:	f2 e9 31 fc ff ff    	bnd jmp 9020 <_init+0x20>
    93ef:	90                   	nop
    93f0:	f3 0f 1e fa          	endbr64 
    93f4:	68 3c 00 00 00       	push   $0x3c
    93f9:	f2 e9 21 fc ff ff    	bnd jmp 9020 <_init+0x20>
    93ff:	90                   	nop
    9400:	f3 0f 1e fa          	endbr64 
    9404:	68 3d 00 00 00       	push   $0x3d
    9409:	f2 e9 11 fc ff ff    	bnd jmp 9020 <_init+0x20>
    940f:	90                   	nop
    9410:	f3 0f 1e fa          	endbr64 
    9414:	68 3e 00 00 00       	push   $0x3e
    9419:	f2 e9 01 fc ff ff    	bnd jmp 9020 <_init+0x20>
    941f:	90                   	nop
    9420:	f3 0f 1e fa          	endbr64 
    9424:	68 3f 00 00 00       	push   $0x3f
    9429:	f2 e9 f1 fb ff ff    	bnd jmp 9020 <_init+0x20>
    942f:	90                   	nop
    9430:	f3 0f 1e fa          	endbr64 
    9434:	68 40 00 00 00       	push   $0x40
    9439:	f2 e9 e1 fb ff ff    	bnd jmp 9020 <_init+0x20>
    943f:	90                   	nop
    9440:	f3 0f 1e fa          	endbr64 
    9444:	68 41 00 00 00       	push   $0x41
    9449:	f2 e9 d1 fb ff ff    	bnd jmp 9020 <_init+0x20>
    944f:	90                   	nop
    9450:	f3 0f 1e fa          	endbr64 
    9454:	68 42 00 00 00       	push   $0x42
    9459:	f2 e9 c1 fb ff ff    	bnd jmp 9020 <_init+0x20>
    945f:	90                   	nop
    9460:	f3 0f 1e fa          	endbr64 
    9464:	68 43 00 00 00       	push   $0x43
    9469:	f2 e9 b1 fb ff ff    	bnd jmp 9020 <_init+0x20>
    946f:	90                   	nop
    9470:	f3 0f 1e fa          	endbr64 
    9474:	68 44 00 00 00       	push   $0x44
    9479:	f2 e9 a1 fb ff ff    	bnd jmp 9020 <_init+0x20>
    947f:	90                   	nop
    9480:	f3 0f 1e fa          	endbr64 
    9484:	68 45 00 00 00       	push   $0x45
    9489:	f2 e9 91 fb ff ff    	bnd jmp 9020 <_init+0x20>
    948f:	90                   	nop
    9490:	f3 0f 1e fa          	endbr64 
    9494:	68 46 00 00 00       	push   $0x46
    9499:	f2 e9 81 fb ff ff    	bnd jmp 9020 <_init+0x20>
    949f:	90                   	nop
    94a0:	f3 0f 1e fa          	endbr64 
    94a4:	68 47 00 00 00       	push   $0x47
    94a9:	f2 e9 71 fb ff ff    	bnd jmp 9020 <_init+0x20>
    94af:	90                   	nop
    94b0:	f3 0f 1e fa          	endbr64 
    94b4:	68 48 00 00 00       	push   $0x48
    94b9:	f2 e9 61 fb ff ff    	bnd jmp 9020 <_init+0x20>
    94bf:	90                   	nop
    94c0:	f3 0f 1e fa          	endbr64 
    94c4:	68 49 00 00 00       	push   $0x49
    94c9:	f2 e9 51 fb ff ff    	bnd jmp 9020 <_init+0x20>
    94cf:	90                   	nop
    94d0:	f3 0f 1e fa          	endbr64 
    94d4:	68 4a 00 00 00       	push   $0x4a
    94d9:	f2 e9 41 fb ff ff    	bnd jmp 9020 <_init+0x20>
    94df:	90                   	nop
    94e0:	f3 0f 1e fa          	endbr64 
    94e4:	68 4b 00 00 00       	push   $0x4b
    94e9:	f2 e9 31 fb ff ff    	bnd jmp 9020 <_init+0x20>
    94ef:	90                   	nop
    94f0:	f3 0f 1e fa          	endbr64 
    94f4:	68 4c 00 00 00       	push   $0x4c
    94f9:	f2 e9 21 fb ff ff    	bnd jmp 9020 <_init+0x20>
    94ff:	90                   	nop
    9500:	f3 0f 1e fa          	endbr64 
    9504:	68 4d 00 00 00       	push   $0x4d
    9509:	f2 e9 11 fb ff ff    	bnd jmp 9020 <_init+0x20>
    950f:	90                   	nop
    9510:	f3 0f 1e fa          	endbr64 
    9514:	68 4e 00 00 00       	push   $0x4e
    9519:	f2 e9 01 fb ff ff    	bnd jmp 9020 <_init+0x20>
    951f:	90                   	nop
    9520:	f3 0f 1e fa          	endbr64 
    9524:	68 4f 00 00 00       	push   $0x4f
    9529:	f2 e9 f1 fa ff ff    	bnd jmp 9020 <_init+0x20>
    952f:	90                   	nop
    9530:	f3 0f 1e fa          	endbr64 
    9534:	68 50 00 00 00       	push   $0x50
    9539:	f2 e9 e1 fa ff ff    	bnd jmp 9020 <_init+0x20>
    953f:	90                   	nop
    9540:	f3 0f 1e fa          	endbr64 
    9544:	68 51 00 00 00       	push   $0x51
    9549:	f2 e9 d1 fa ff ff    	bnd jmp 9020 <_init+0x20>
    954f:	90                   	nop
    9550:	f3 0f 1e fa          	endbr64 
    9554:	68 52 00 00 00       	push   $0x52
    9559:	f2 e9 c1 fa ff ff    	bnd jmp 9020 <_init+0x20>
    955f:	90                   	nop
    9560:	f3 0f 1e fa          	endbr64 
    9564:	68 53 00 00 00       	push   $0x53
    9569:	f2 e9 b1 fa ff ff    	bnd jmp 9020 <_init+0x20>
    956f:	90                   	nop
    9570:	f3 0f 1e fa          	endbr64 
    9574:	68 54 00 00 00       	push   $0x54
    9579:	f2 e9 a1 fa ff ff    	bnd jmp 9020 <_init+0x20>
    957f:	90                   	nop
    9580:	f3 0f 1e fa          	endbr64 
    9584:	68 55 00 00 00       	push   $0x55
    9589:	f2 e9 91 fa ff ff    	bnd jmp 9020 <_init+0x20>
    958f:	90                   	nop
    9590:	f3 0f 1e fa          	endbr64 
    9594:	68 56 00 00 00       	push   $0x56
    9599:	f2 e9 81 fa ff ff    	bnd jmp 9020 <_init+0x20>
    959f:	90                   	nop
    95a0:	f3 0f 1e fa          	endbr64 
    95a4:	68 57 00 00 00       	push   $0x57
    95a9:	f2 e9 71 fa ff ff    	bnd jmp 9020 <_init+0x20>
    95af:	90                   	nop
    95b0:	f3 0f 1e fa          	endbr64 
    95b4:	68 58 00 00 00       	push   $0x58
    95b9:	f2 e9 61 fa ff ff    	bnd jmp 9020 <_init+0x20>
    95bf:	90                   	nop
    95c0:	f3 0f 1e fa          	endbr64 
    95c4:	68 59 00 00 00       	push   $0x59
    95c9:	f2 e9 51 fa ff ff    	bnd jmp 9020 <_init+0x20>
    95cf:	90                   	nop
    95d0:	f3 0f 1e fa          	endbr64 
    95d4:	68 5a 00 00 00       	push   $0x5a
    95d9:	f2 e9 41 fa ff ff    	bnd jmp 9020 <_init+0x20>
    95df:	90                   	nop
    95e0:	f3 0f 1e fa          	endbr64 
    95e4:	68 5b 00 00 00       	push   $0x5b
    95e9:	f2 e9 31 fa ff ff    	bnd jmp 9020 <_init+0x20>
    95ef:	90                   	nop
    95f0:	f3 0f 1e fa          	endbr64 
    95f4:	68 5c 00 00 00       	push   $0x5c
    95f9:	f2 e9 21 fa ff ff    	bnd jmp 9020 <_init+0x20>
    95ff:	90                   	nop
    9600:	f3 0f 1e fa          	endbr64 
    9604:	68 5d 00 00 00       	push   $0x5d
    9609:	f2 e9 11 fa ff ff    	bnd jmp 9020 <_init+0x20>
    960f:	90                   	nop
    9610:	f3 0f 1e fa          	endbr64 
    9614:	68 5e 00 00 00       	push   $0x5e
    9619:	f2 e9 01 fa ff ff    	bnd jmp 9020 <_init+0x20>
    961f:	90                   	nop
    9620:	f3 0f 1e fa          	endbr64 
    9624:	68 5f 00 00 00       	push   $0x5f
    9629:	f2 e9 f1 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    962f:	90                   	nop
    9630:	f3 0f 1e fa          	endbr64 
    9634:	68 60 00 00 00       	push   $0x60
    9639:	f2 e9 e1 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    963f:	90                   	nop
    9640:	f3 0f 1e fa          	endbr64 
    9644:	68 61 00 00 00       	push   $0x61
    9649:	f2 e9 d1 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    964f:	90                   	nop
    9650:	f3 0f 1e fa          	endbr64 
    9654:	68 62 00 00 00       	push   $0x62
    9659:	f2 e9 c1 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    965f:	90                   	nop
    9660:	f3 0f 1e fa          	endbr64 
    9664:	68 63 00 00 00       	push   $0x63
    9669:	f2 e9 b1 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    966f:	90                   	nop
    9670:	f3 0f 1e fa          	endbr64 
    9674:	68 64 00 00 00       	push   $0x64
    9679:	f2 e9 a1 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    967f:	90                   	nop
    9680:	f3 0f 1e fa          	endbr64 
    9684:	68 65 00 00 00       	push   $0x65
    9689:	f2 e9 91 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    968f:	90                   	nop
    9690:	f3 0f 1e fa          	endbr64 
    9694:	68 66 00 00 00       	push   $0x66
    9699:	f2 e9 81 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    969f:	90                   	nop
    96a0:	f3 0f 1e fa          	endbr64 
    96a4:	68 67 00 00 00       	push   $0x67
    96a9:	f2 e9 71 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    96af:	90                   	nop
    96b0:	f3 0f 1e fa          	endbr64 
    96b4:	68 68 00 00 00       	push   $0x68
    96b9:	f2 e9 61 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    96bf:	90                   	nop
    96c0:	f3 0f 1e fa          	endbr64 
    96c4:	68 69 00 00 00       	push   $0x69
    96c9:	f2 e9 51 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    96cf:	90                   	nop
    96d0:	f3 0f 1e fa          	endbr64 
    96d4:	68 6a 00 00 00       	push   $0x6a
    96d9:	f2 e9 41 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    96df:	90                   	nop
    96e0:	f3 0f 1e fa          	endbr64 
    96e4:	68 6b 00 00 00       	push   $0x6b
    96e9:	f2 e9 31 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    96ef:	90                   	nop
    96f0:	f3 0f 1e fa          	endbr64 
    96f4:	68 6c 00 00 00       	push   $0x6c
    96f9:	f2 e9 21 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    96ff:	90                   	nop
    9700:	f3 0f 1e fa          	endbr64 
    9704:	68 6d 00 00 00       	push   $0x6d
    9709:	f2 e9 11 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    970f:	90                   	nop
    9710:	f3 0f 1e fa          	endbr64 
    9714:	68 6e 00 00 00       	push   $0x6e
    9719:	f2 e9 01 f9 ff ff    	bnd jmp 9020 <_init+0x20>
    971f:	90                   	nop
    9720:	f3 0f 1e fa          	endbr64 
    9724:	68 6f 00 00 00       	push   $0x6f
    9729:	f2 e9 f1 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    972f:	90                   	nop
    9730:	f3 0f 1e fa          	endbr64 
    9734:	68 70 00 00 00       	push   $0x70
    9739:	f2 e9 e1 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    973f:	90                   	nop
    9740:	f3 0f 1e fa          	endbr64 
    9744:	68 71 00 00 00       	push   $0x71
    9749:	f2 e9 d1 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    974f:	90                   	nop
    9750:	f3 0f 1e fa          	endbr64 
    9754:	68 72 00 00 00       	push   $0x72
    9759:	f2 e9 c1 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    975f:	90                   	nop
    9760:	f3 0f 1e fa          	endbr64 
    9764:	68 73 00 00 00       	push   $0x73
    9769:	f2 e9 b1 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    976f:	90                   	nop
    9770:	f3 0f 1e fa          	endbr64 
    9774:	68 74 00 00 00       	push   $0x74
    9779:	f2 e9 a1 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    977f:	90                   	nop
    9780:	f3 0f 1e fa          	endbr64 
    9784:	68 75 00 00 00       	push   $0x75
    9789:	f2 e9 91 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    978f:	90                   	nop
    9790:	f3 0f 1e fa          	endbr64 
    9794:	68 76 00 00 00       	push   $0x76
    9799:	f2 e9 81 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    979f:	90                   	nop
    97a0:	f3 0f 1e fa          	endbr64 
    97a4:	68 77 00 00 00       	push   $0x77
    97a9:	f2 e9 71 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    97af:	90                   	nop
    97b0:	f3 0f 1e fa          	endbr64 
    97b4:	68 78 00 00 00       	push   $0x78
    97b9:	f2 e9 61 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    97bf:	90                   	nop
    97c0:	f3 0f 1e fa          	endbr64 
    97c4:	68 79 00 00 00       	push   $0x79
    97c9:	f2 e9 51 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    97cf:	90                   	nop
    97d0:	f3 0f 1e fa          	endbr64 
    97d4:	68 7a 00 00 00       	push   $0x7a
    97d9:	f2 e9 41 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    97df:	90                   	nop
    97e0:	f3 0f 1e fa          	endbr64 
    97e4:	68 7b 00 00 00       	push   $0x7b
    97e9:	f2 e9 31 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    97ef:	90                   	nop
    97f0:	f3 0f 1e fa          	endbr64 
    97f4:	68 7c 00 00 00       	push   $0x7c
    97f9:	f2 e9 21 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    97ff:	90                   	nop
    9800:	f3 0f 1e fa          	endbr64 
    9804:	68 7d 00 00 00       	push   $0x7d
    9809:	f2 e9 11 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    980f:	90                   	nop
    9810:	f3 0f 1e fa          	endbr64 
    9814:	68 7e 00 00 00       	push   $0x7e
    9819:	f2 e9 01 f8 ff ff    	bnd jmp 9020 <_init+0x20>
    981f:	90                   	nop
    9820:	f3 0f 1e fa          	endbr64 
    9824:	68 7f 00 00 00       	push   $0x7f
    9829:	f2 e9 f1 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    982f:	90                   	nop
    9830:	f3 0f 1e fa          	endbr64 
    9834:	68 80 00 00 00       	push   $0x80
    9839:	f2 e9 e1 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    983f:	90                   	nop
    9840:	f3 0f 1e fa          	endbr64 
    9844:	68 81 00 00 00       	push   $0x81
    9849:	f2 e9 d1 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    984f:	90                   	nop
    9850:	f3 0f 1e fa          	endbr64 
    9854:	68 82 00 00 00       	push   $0x82
    9859:	f2 e9 c1 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    985f:	90                   	nop
    9860:	f3 0f 1e fa          	endbr64 
    9864:	68 83 00 00 00       	push   $0x83
    9869:	f2 e9 b1 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    986f:	90                   	nop
    9870:	f3 0f 1e fa          	endbr64 
    9874:	68 84 00 00 00       	push   $0x84
    9879:	f2 e9 a1 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    987f:	90                   	nop
    9880:	f3 0f 1e fa          	endbr64 
    9884:	68 85 00 00 00       	push   $0x85
    9889:	f2 e9 91 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    988f:	90                   	nop
    9890:	f3 0f 1e fa          	endbr64 
    9894:	68 86 00 00 00       	push   $0x86
    9899:	f2 e9 81 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    989f:	90                   	nop
    98a0:	f3 0f 1e fa          	endbr64 
    98a4:	68 87 00 00 00       	push   $0x87
    98a9:	f2 e9 71 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    98af:	90                   	nop
    98b0:	f3 0f 1e fa          	endbr64 
    98b4:	68 88 00 00 00       	push   $0x88
    98b9:	f2 e9 61 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    98bf:	90                   	nop
    98c0:	f3 0f 1e fa          	endbr64 
    98c4:	68 89 00 00 00       	push   $0x89
    98c9:	f2 e9 51 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    98cf:	90                   	nop
    98d0:	f3 0f 1e fa          	endbr64 
    98d4:	68 8a 00 00 00       	push   $0x8a
    98d9:	f2 e9 41 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    98df:	90                   	nop
    98e0:	f3 0f 1e fa          	endbr64 
    98e4:	68 8b 00 00 00       	push   $0x8b
    98e9:	f2 e9 31 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    98ef:	90                   	nop
    98f0:	f3 0f 1e fa          	endbr64 
    98f4:	68 8c 00 00 00       	push   $0x8c
    98f9:	f2 e9 21 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    98ff:	90                   	nop
    9900:	f3 0f 1e fa          	endbr64 
    9904:	68 8d 00 00 00       	push   $0x8d
    9909:	f2 e9 11 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    990f:	90                   	nop
    9910:	f3 0f 1e fa          	endbr64 
    9914:	68 8e 00 00 00       	push   $0x8e
    9919:	f2 e9 01 f7 ff ff    	bnd jmp 9020 <_init+0x20>
    991f:	90                   	nop
    9920:	f3 0f 1e fa          	endbr64 
    9924:	68 8f 00 00 00       	push   $0x8f
    9929:	f2 e9 f1 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    992f:	90                   	nop
    9930:	f3 0f 1e fa          	endbr64 
    9934:	68 90 00 00 00       	push   $0x90
    9939:	f2 e9 e1 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    993f:	90                   	nop
    9940:	f3 0f 1e fa          	endbr64 
    9944:	68 91 00 00 00       	push   $0x91
    9949:	f2 e9 d1 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    994f:	90                   	nop
    9950:	f3 0f 1e fa          	endbr64 
    9954:	68 92 00 00 00       	push   $0x92
    9959:	f2 e9 c1 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    995f:	90                   	nop
    9960:	f3 0f 1e fa          	endbr64 
    9964:	68 93 00 00 00       	push   $0x93
    9969:	f2 e9 b1 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    996f:	90                   	nop
    9970:	f3 0f 1e fa          	endbr64 
    9974:	68 94 00 00 00       	push   $0x94
    9979:	f2 e9 a1 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    997f:	90                   	nop
    9980:	f3 0f 1e fa          	endbr64 
    9984:	68 95 00 00 00       	push   $0x95
    9989:	f2 e9 91 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    998f:	90                   	nop
    9990:	f3 0f 1e fa          	endbr64 
    9994:	68 96 00 00 00       	push   $0x96
    9999:	f2 e9 81 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    999f:	90                   	nop
    99a0:	f3 0f 1e fa          	endbr64 
    99a4:	68 97 00 00 00       	push   $0x97
    99a9:	f2 e9 71 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    99af:	90                   	nop
    99b0:	f3 0f 1e fa          	endbr64 
    99b4:	68 98 00 00 00       	push   $0x98
    99b9:	f2 e9 61 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    99bf:	90                   	nop
    99c0:	f3 0f 1e fa          	endbr64 
    99c4:	68 99 00 00 00       	push   $0x99
    99c9:	f2 e9 51 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    99cf:	90                   	nop
    99d0:	f3 0f 1e fa          	endbr64 
    99d4:	68 9a 00 00 00       	push   $0x9a
    99d9:	f2 e9 41 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    99df:	90                   	nop
    99e0:	f3 0f 1e fa          	endbr64 
    99e4:	68 9b 00 00 00       	push   $0x9b
    99e9:	f2 e9 31 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    99ef:	90                   	nop
    99f0:	f3 0f 1e fa          	endbr64 
    99f4:	68 9c 00 00 00       	push   $0x9c
    99f9:	f2 e9 21 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    99ff:	90                   	nop
    9a00:	f3 0f 1e fa          	endbr64 
    9a04:	68 9d 00 00 00       	push   $0x9d
    9a09:	f2 e9 11 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    9a0f:	90                   	nop
    9a10:	f3 0f 1e fa          	endbr64 
    9a14:	68 9e 00 00 00       	push   $0x9e
    9a19:	f2 e9 01 f6 ff ff    	bnd jmp 9020 <_init+0x20>
    9a1f:	90                   	nop
    9a20:	f3 0f 1e fa          	endbr64 
    9a24:	68 9f 00 00 00       	push   $0x9f
    9a29:	f2 e9 f1 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9a2f:	90                   	nop
    9a30:	f3 0f 1e fa          	endbr64 
    9a34:	68 a0 00 00 00       	push   $0xa0
    9a39:	f2 e9 e1 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9a3f:	90                   	nop
    9a40:	f3 0f 1e fa          	endbr64 
    9a44:	68 a1 00 00 00       	push   $0xa1
    9a49:	f2 e9 d1 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9a4f:	90                   	nop
    9a50:	f3 0f 1e fa          	endbr64 
    9a54:	68 a2 00 00 00       	push   $0xa2
    9a59:	f2 e9 c1 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9a5f:	90                   	nop
    9a60:	f3 0f 1e fa          	endbr64 
    9a64:	68 a3 00 00 00       	push   $0xa3
    9a69:	f2 e9 b1 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9a6f:	90                   	nop
    9a70:	f3 0f 1e fa          	endbr64 
    9a74:	68 a4 00 00 00       	push   $0xa4
    9a79:	f2 e9 a1 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9a7f:	90                   	nop
    9a80:	f3 0f 1e fa          	endbr64 
    9a84:	68 a5 00 00 00       	push   $0xa5
    9a89:	f2 e9 91 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9a8f:	90                   	nop
    9a90:	f3 0f 1e fa          	endbr64 
    9a94:	68 a6 00 00 00       	push   $0xa6
    9a99:	f2 e9 81 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9a9f:	90                   	nop
    9aa0:	f3 0f 1e fa          	endbr64 
    9aa4:	68 a7 00 00 00       	push   $0xa7
    9aa9:	f2 e9 71 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9aaf:	90                   	nop
    9ab0:	f3 0f 1e fa          	endbr64 
    9ab4:	68 a8 00 00 00       	push   $0xa8
    9ab9:	f2 e9 61 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9abf:	90                   	nop
    9ac0:	f3 0f 1e fa          	endbr64 
    9ac4:	68 a9 00 00 00       	push   $0xa9
    9ac9:	f2 e9 51 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9acf:	90                   	nop
    9ad0:	f3 0f 1e fa          	endbr64 
    9ad4:	68 aa 00 00 00       	push   $0xaa
    9ad9:	f2 e9 41 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9adf:	90                   	nop
    9ae0:	f3 0f 1e fa          	endbr64 
    9ae4:	68 ab 00 00 00       	push   $0xab
    9ae9:	f2 e9 31 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9aef:	90                   	nop
    9af0:	f3 0f 1e fa          	endbr64 
    9af4:	68 ac 00 00 00       	push   $0xac
    9af9:	f2 e9 21 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9aff:	90                   	nop
    9b00:	f3 0f 1e fa          	endbr64 
    9b04:	68 ad 00 00 00       	push   $0xad
    9b09:	f2 e9 11 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9b0f:	90                   	nop
    9b10:	f3 0f 1e fa          	endbr64 
    9b14:	68 ae 00 00 00       	push   $0xae
    9b19:	f2 e9 01 f5 ff ff    	bnd jmp 9020 <_init+0x20>
    9b1f:	90                   	nop
    9b20:	f3 0f 1e fa          	endbr64 
    9b24:	68 af 00 00 00       	push   $0xaf
    9b29:	f2 e9 f1 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9b2f:	90                   	nop
    9b30:	f3 0f 1e fa          	endbr64 
    9b34:	68 b0 00 00 00       	push   $0xb0
    9b39:	f2 e9 e1 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9b3f:	90                   	nop
    9b40:	f3 0f 1e fa          	endbr64 
    9b44:	68 b1 00 00 00       	push   $0xb1
    9b49:	f2 e9 d1 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9b4f:	90                   	nop
    9b50:	f3 0f 1e fa          	endbr64 
    9b54:	68 b2 00 00 00       	push   $0xb2
    9b59:	f2 e9 c1 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9b5f:	90                   	nop
    9b60:	f3 0f 1e fa          	endbr64 
    9b64:	68 b3 00 00 00       	push   $0xb3
    9b69:	f2 e9 b1 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9b6f:	90                   	nop
    9b70:	f3 0f 1e fa          	endbr64 
    9b74:	68 b4 00 00 00       	push   $0xb4
    9b79:	f2 e9 a1 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9b7f:	90                   	nop
    9b80:	f3 0f 1e fa          	endbr64 
    9b84:	68 b5 00 00 00       	push   $0xb5
    9b89:	f2 e9 91 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9b8f:	90                   	nop
    9b90:	f3 0f 1e fa          	endbr64 
    9b94:	68 b6 00 00 00       	push   $0xb6
    9b99:	f2 e9 81 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9b9f:	90                   	nop
    9ba0:	f3 0f 1e fa          	endbr64 
    9ba4:	68 b7 00 00 00       	push   $0xb7
    9ba9:	f2 e9 71 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9baf:	90                   	nop
    9bb0:	f3 0f 1e fa          	endbr64 
    9bb4:	68 b8 00 00 00       	push   $0xb8
    9bb9:	f2 e9 61 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9bbf:	90                   	nop
    9bc0:	f3 0f 1e fa          	endbr64 
    9bc4:	68 b9 00 00 00       	push   $0xb9
    9bc9:	f2 e9 51 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9bcf:	90                   	nop
    9bd0:	f3 0f 1e fa          	endbr64 
    9bd4:	68 ba 00 00 00       	push   $0xba
    9bd9:	f2 e9 41 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9bdf:	90                   	nop
    9be0:	f3 0f 1e fa          	endbr64 
    9be4:	68 bb 00 00 00       	push   $0xbb
    9be9:	f2 e9 31 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9bef:	90                   	nop
    9bf0:	f3 0f 1e fa          	endbr64 
    9bf4:	68 bc 00 00 00       	push   $0xbc
    9bf9:	f2 e9 21 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9bff:	90                   	nop
    9c00:	f3 0f 1e fa          	endbr64 
    9c04:	68 bd 00 00 00       	push   $0xbd
    9c09:	f2 e9 11 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9c0f:	90                   	nop
    9c10:	f3 0f 1e fa          	endbr64 
    9c14:	68 be 00 00 00       	push   $0xbe
    9c19:	f2 e9 01 f4 ff ff    	bnd jmp 9020 <_init+0x20>
    9c1f:	90                   	nop
    9c20:	f3 0f 1e fa          	endbr64 
    9c24:	68 bf 00 00 00       	push   $0xbf
    9c29:	f2 e9 f1 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9c2f:	90                   	nop
    9c30:	f3 0f 1e fa          	endbr64 
    9c34:	68 c0 00 00 00       	push   $0xc0
    9c39:	f2 e9 e1 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9c3f:	90                   	nop
    9c40:	f3 0f 1e fa          	endbr64 
    9c44:	68 c1 00 00 00       	push   $0xc1
    9c49:	f2 e9 d1 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9c4f:	90                   	nop
    9c50:	f3 0f 1e fa          	endbr64 
    9c54:	68 c2 00 00 00       	push   $0xc2
    9c59:	f2 e9 c1 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9c5f:	90                   	nop
    9c60:	f3 0f 1e fa          	endbr64 
    9c64:	68 c3 00 00 00       	push   $0xc3
    9c69:	f2 e9 b1 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9c6f:	90                   	nop
    9c70:	f3 0f 1e fa          	endbr64 
    9c74:	68 c4 00 00 00       	push   $0xc4
    9c79:	f2 e9 a1 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9c7f:	90                   	nop
    9c80:	f3 0f 1e fa          	endbr64 
    9c84:	68 c5 00 00 00       	push   $0xc5
    9c89:	f2 e9 91 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9c8f:	90                   	nop
    9c90:	f3 0f 1e fa          	endbr64 
    9c94:	68 c6 00 00 00       	push   $0xc6
    9c99:	f2 e9 81 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9c9f:	90                   	nop
    9ca0:	f3 0f 1e fa          	endbr64 
    9ca4:	68 c7 00 00 00       	push   $0xc7
    9ca9:	f2 e9 71 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9caf:	90                   	nop
    9cb0:	f3 0f 1e fa          	endbr64 
    9cb4:	68 c8 00 00 00       	push   $0xc8
    9cb9:	f2 e9 61 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9cbf:	90                   	nop
    9cc0:	f3 0f 1e fa          	endbr64 
    9cc4:	68 c9 00 00 00       	push   $0xc9
    9cc9:	f2 e9 51 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9ccf:	90                   	nop
    9cd0:	f3 0f 1e fa          	endbr64 
    9cd4:	68 ca 00 00 00       	push   $0xca
    9cd9:	f2 e9 41 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9cdf:	90                   	nop
    9ce0:	f3 0f 1e fa          	endbr64 
    9ce4:	68 cb 00 00 00       	push   $0xcb
    9ce9:	f2 e9 31 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9cef:	90                   	nop
    9cf0:	f3 0f 1e fa          	endbr64 
    9cf4:	68 cc 00 00 00       	push   $0xcc
    9cf9:	f2 e9 21 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9cff:	90                   	nop
    9d00:	f3 0f 1e fa          	endbr64 
    9d04:	68 cd 00 00 00       	push   $0xcd
    9d09:	f2 e9 11 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9d0f:	90                   	nop
    9d10:	f3 0f 1e fa          	endbr64 
    9d14:	68 ce 00 00 00       	push   $0xce
    9d19:	f2 e9 01 f3 ff ff    	bnd jmp 9020 <_init+0x20>
    9d1f:	90                   	nop
    9d20:	f3 0f 1e fa          	endbr64 
    9d24:	68 cf 00 00 00       	push   $0xcf
    9d29:	f2 e9 f1 f2 ff ff    	bnd jmp 9020 <_init+0x20>
    9d2f:	90                   	nop
    9d30:	f3 0f 1e fa          	endbr64 
    9d34:	68 d0 00 00 00       	push   $0xd0
    9d39:	f2 e9 e1 f2 ff ff    	bnd jmp 9020 <_init+0x20>
    9d3f:	90                   	nop
    9d40:	f3 0f 1e fa          	endbr64 
    9d44:	68 d1 00 00 00       	push   $0xd1
    9d49:	f2 e9 d1 f2 ff ff    	bnd jmp 9020 <_init+0x20>
    9d4f:	90                   	nop
    9d50:	f3 0f 1e fa          	endbr64 
    9d54:	68 d2 00 00 00       	push   $0xd2
    9d59:	f2 e9 c1 f2 ff ff    	bnd jmp 9020 <_init+0x20>
    9d5f:	90                   	nop
    9d60:	f3 0f 1e fa          	endbr64 
    9d64:	68 d3 00 00 00       	push   $0xd3
    9d69:	f2 e9 b1 f2 ff ff    	bnd jmp 9020 <_init+0x20>
    9d6f:	90                   	nop
    9d70:	f3 0f 1e fa          	endbr64 
    9d74:	68 d4 00 00 00       	push   $0xd4
    9d79:	f2 e9 a1 f2 ff ff    	bnd jmp 9020 <_init+0x20>
    9d7f:	90                   	nop
    9d80:	f3 0f 1e fa          	endbr64 
    9d84:	68 d5 00 00 00       	push   $0xd5
    9d89:	f2 e9 91 f2 ff ff    	bnd jmp 9020 <_init+0x20>
    9d8f:	90                   	nop
    9d90:	f3 0f 1e fa          	endbr64 
    9d94:	68 d6 00 00 00       	push   $0xd6
    9d99:	f2 e9 81 f2 ff ff    	bnd jmp 9020 <_init+0x20>
    9d9f:	90                   	nop
    9da0:	f3 0f 1e fa          	endbr64 
    9da4:	68 d7 00 00 00       	push   $0xd7
    9da9:	f2 e9 71 f2 ff ff    	bnd jmp 9020 <_init+0x20>
    9daf:	90                   	nop

Disassembly of section .plt.got:

0000000000009db0 <__cxa_finalize@plt>:
    9db0:	f3 0f 1e fa          	endbr64 
    9db4:	f2 ff 25 0d 32 01 00 	bnd jmp *0x1320d(%rip)        # 1cfc8 <__cxa_finalize@GLIBC_2.2.5>
    9dbb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

Disassembly of section .plt.sec:

0000000000009dc0 <_Znam@plt>:
    9dc0:	f3 0f 1e fa          	endbr64 
    9dc4:	f2 ff 25 35 2b 01 00 	bnd jmp *0x12b35(%rip)        # 1c900 <_Znam@GLIBCXX_3.4>
    9dcb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009dd0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE8capacityEv@plt>:
    9dd0:	f3 0f 1e fa          	endbr64 
    9dd4:	f2 ff 25 2d 2b 01 00 	bnd jmp *0x12b2d(%rip)        # 1c908 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE8capacityEv@GLIBCXX_3.4.21>
    9ddb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009de0 <printf@plt>:
    9de0:	f3 0f 1e fa          	endbr64 
    9de4:	f2 ff 25 25 2b 01 00 	bnd jmp *0x12b25(%rip)        # 1c910 <printf@GLIBC_2.2.5>
    9deb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009df0 <MPI_Win_unlock@plt>:
    9df0:	f3 0f 1e fa          	endbr64 
    9df4:	f2 ff 25 1d 2b 01 00 	bnd jmp *0x12b1d(%rip)        # 1c918 <MPI_Win_unlock@Base>
    9dfb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e00 <MPI_Group_excl@plt>:
    9e00:	f3 0f 1e fa          	endbr64 
    9e04:	f2 ff 25 15 2b 01 00 	bnd jmp *0x12b15(%rip)        # 1c920 <MPI_Group_excl@Base>
    9e0b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e10 <MPI_Allgather@plt>:
    9e10:	f3 0f 1e fa          	endbr64 
    9e14:	f2 ff 25 0d 2b 01 00 	bnd jmp *0x12b0d(%rip)        # 1c928 <MPI_Allgather@Base>
    9e1b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e20 <MPI_Recv_init@plt>:
    9e20:	f3 0f 1e fa          	endbr64 
    9e24:	f2 ff 25 05 2b 01 00 	bnd jmp *0x12b05(%rip)        # 1c930 <MPI_Recv_init@Base>
    9e2b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e30 <MPI_Cart_create@plt>:
    9e30:	f3 0f 1e fa          	endbr64 
    9e34:	f2 ff 25 fd 2a 01 00 	bnd jmp *0x12afd(%rip)        # 1c938 <MPI_Cart_create@Base>
    9e3b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e40 <MPI_Comm_get_name@plt>:
    9e40:	f3 0f 1e fa          	endbr64 
    9e44:	f2 ff 25 f5 2a 01 00 	bnd jmp *0x12af5(%rip)        # 1c940 <MPI_Comm_get_name@Base>
    9e4b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e50 <MPI_Info_dup@plt>:
    9e50:	f3 0f 1e fa          	endbr64 
    9e54:	f2 ff 25 ed 2a 01 00 	bnd jmp *0x12aed(%rip)        # 1c948 <MPI_Info_dup@Base>
    9e5b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e60 <MPI_Irsend@plt>:
    9e60:	f3 0f 1e fa          	endbr64 
    9e64:	f2 ff 25 e5 2a 01 00 	bnd jmp *0x12ae5(%rip)        # 1c950 <MPI_Irsend@Base>
    9e6b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e70 <MPI_Bcast@plt>:
    9e70:	f3 0f 1e fa          	endbr64 
    9e74:	f2 ff 25 dd 2a 01 00 	bnd jmp *0x12add(%rip)        # 1c958 <MPI_Bcast@Base>
    9e7b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e80 <_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE@plt>:
    9e80:	f3 0f 1e fa          	endbr64 
    9e84:	f2 ff 25 d5 2a 01 00 	bnd jmp *0x12ad5(%rip)        # 1c960 <_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE@GLIBCXX_3.4.21>
    9e8b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e90 <MPI_Ssend@plt>:
    9e90:	f3 0f 1e fa          	endbr64 
    9e94:	f2 ff 25 cd 2a 01 00 	bnd jmp *0x12acd(%rip)        # 1c968 <MPI_Ssend@Base>
    9e9b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009ea0 <MPI_Group_range_incl@plt>:
    9ea0:	f3 0f 1e fa          	endbr64 
    9ea4:	f2 ff 25 c5 2a 01 00 	bnd jmp *0x12ac5(%rip)        # 1c970 <MPI_Group_range_incl@Base>
    9eab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009eb0 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEE5closeEv@plt>:
    9eb0:	f3 0f 1e fa          	endbr64 
    9eb4:	f2 ff 25 bd 2a 01 00 	bnd jmp *0x12abd(%rip)        # 1c978 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEE5closeEv@GLIBCXX_3.4>
    9ebb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009ec0 <MPI_Op_create@plt>:
    9ec0:	f3 0f 1e fa          	endbr64 
    9ec4:	f2 ff 25 b5 2a 01 00 	bnd jmp *0x12ab5(%rip)        # 1c980 <MPI_Op_create@Base>
    9ecb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009ed0 <MPI_Group_size@plt>:
    9ed0:	f3 0f 1e fa          	endbr64 
    9ed4:	f2 ff 25 ad 2a 01 00 	bnd jmp *0x12aad(%rip)        # 1c988 <MPI_Group_size@Base>
    9edb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009ee0 <MPI_Comm_disconnect@plt>:
    9ee0:	f3 0f 1e fa          	endbr64 
    9ee4:	f2 ff 25 a5 2a 01 00 	bnd jmp *0x12aa5(%rip)        # 1c990 <MPI_Comm_disconnect@Base>
    9eeb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009ef0 <_ZNSolsEm@plt>:
    9ef0:	f3 0f 1e fa          	endbr64 
    9ef4:	f2 ff 25 9d 2a 01 00 	bnd jmp *0x12a9d(%rip)        # 1c998 <_ZNSolsEm@GLIBCXX_3.4>
    9efb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f00 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmPKc@plt>:
    9f00:	f3 0f 1e fa          	endbr64 
    9f04:	f2 ff 25 95 2a 01 00 	bnd jmp *0x12a95(%rip)        # 1c9a0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmPKc@GLIBCXX_3.4.21>
    9f0b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f10 <MPI_Win_get_name@plt>:
    9f10:	f3 0f 1e fa          	endbr64 
    9f14:	f2 ff 25 8d 2a 01 00 	bnd jmp *0x12a8d(%rip)        # 1c9a8 <MPI_Win_get_name@Base>
    9f1b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f20 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEED1Ev@plt>:
    9f20:	f3 0f 1e fa          	endbr64 
    9f24:	f2 ff 25 85 2a 01 00 	bnd jmp *0x12a85(%rip)        # 1c9b0 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEED1Ev@GLIBCXX_3.4>
    9f2b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f30 <MPI_Pack_size@plt>:
    9f30:	f3 0f 1e fa          	endbr64 
    9f34:	f2 ff 25 7d 2a 01 00 	bnd jmp *0x12a7d(%rip)        # 1c9b8 <MPI_Pack_size@Base>
    9f3b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f40 <MPI_Request_get_status@plt>:
    9f40:	f3 0f 1e fa          	endbr64 
    9f44:	f2 ff 25 75 2a 01 00 	bnd jmp *0x12a75(%rip)        # 1c9c0 <MPI_Request_get_status@Base>
    9f4b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f50 <MPI_Graph_neighbors@plt>:
    9f50:	f3 0f 1e fa          	endbr64 
    9f54:	f2 ff 25 6d 2a 01 00 	bnd jmp *0x12a6d(%rip)        # 1c9c8 <MPI_Graph_neighbors@Base>
    9f5b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f60 <MPI_Type_create_hvector@plt>:
    9f60:	f3 0f 1e fa          	endbr64 
    9f64:	f2 ff 25 65 2a 01 00 	bnd jmp *0x12a65(%rip)        # 1c9d0 <MPI_Type_create_hvector@Base>
    9f6b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f70 <MPI_Unpack_external@plt>:
    9f70:	f3 0f 1e fa          	endbr64 
    9f74:	f2 ff 25 5d 2a 01 00 	bnd jmp *0x12a5d(%rip)        # 1c9d8 <MPI_Unpack_external@Base>
    9f7b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f80 <MPI_Probe@plt>:
    9f80:	f3 0f 1e fa          	endbr64 
    9f84:	f2 ff 25 55 2a 01 00 	bnd jmp *0x12a55(%rip)        # 1c9e0 <MPI_Probe@Base>
    9f8b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f90 <_ZNSolsEPFRSt8ios_baseS0_E@plt>:
    9f90:	f3 0f 1e fa          	endbr64 
    9f94:	f2 ff 25 4d 2a 01 00 	bnd jmp *0x12a4d(%rip)        # 1c9e8 <_ZNSolsEPFRSt8ios_baseS0_E@GLIBCXX_3.4>
    9f9b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009fa0 <MPI_Send@plt>:
    9fa0:	f3 0f 1e fa          	endbr64 
    9fa4:	f2 ff 25 45 2a 01 00 	bnd jmp *0x12a45(%rip)        # 1c9f0 <MPI_Send@Base>
    9fab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009fb0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE13get_allocatorEv@plt>:
    9fb0:	f3 0f 1e fa          	endbr64 
    9fb4:	f2 ff 25 3d 2a 01 00 	bnd jmp *0x12a3d(%rip)        # 1c9f8 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE13get_allocatorEv@GLIBCXX_3.4.21>
    9fbb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009fc0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE5c_strEv@plt>:
    9fc0:	f3 0f 1e fa          	endbr64 
    9fc4:	f2 ff 25 35 2a 01 00 	bnd jmp *0x12a35(%rip)        # 1ca00 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE5c_strEv@GLIBCXX_3.4.21>
    9fcb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009fd0 <MPI_Info_set@plt>:
    9fd0:	f3 0f 1e fa          	endbr64 
    9fd4:	f2 ff 25 2d 2a 01 00 	bnd jmp *0x12a2d(%rip)        # 1ca08 <MPI_Info_set@Base>
    9fdb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009fe0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>:
    9fe0:	f3 0f 1e fa          	endbr64 
    9fe4:	f2 ff 25 25 2a 01 00 	bnd jmp *0x12a25(%rip)        # 1ca10 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@GLIBCXX_3.4.21>
    9feb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009ff0 <MPI_Type_get_extent@plt>:
    9ff0:	f3 0f 1e fa          	endbr64 
    9ff4:	f2 ff 25 1d 2a 01 00 	bnd jmp *0x12a1d(%rip)        # 1ca18 <MPI_Type_get_extent@Base>
    9ffb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a000 <MPI_Win_post@plt>:
    a000:	f3 0f 1e fa          	endbr64 
    a004:	f2 ff 25 15 2a 01 00 	bnd jmp *0x12a15(%rip)        # 1ca20 <MPI_Win_post@Base>
    a00b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a010 <MPI_Isend@plt>:
    a010:	f3 0f 1e fa          	endbr64 
    a014:	f2 ff 25 0d 2a 01 00 	bnd jmp *0x12a0d(%rip)        # 1ca28 <MPI_Isend@Base>
    a01b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a020 <MPI_Unpack@plt>:
    a020:	f3 0f 1e fa          	endbr64 
    a024:	f2 ff 25 05 2a 01 00 	bnd jmp *0x12a05(%rip)        # 1ca30 <MPI_Unpack@Base>
    a02b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a030 <MPI_Iprobe@plt>:
    a030:	f3 0f 1e fa          	endbr64 
    a034:	f2 ff 25 fd 29 01 00 	bnd jmp *0x129fd(%rip)        # 1ca38 <MPI_Iprobe@Base>
    a03b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a040 <MPI_Comm_set_attr@plt>:
    a040:	f3 0f 1e fa          	endbr64 
    a044:	f2 ff 25 f5 29 01 00 	bnd jmp *0x129f5(%rip)        # 1ca40 <MPI_Comm_set_attr@Base>
    a04b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a050 <MPI_Type_get_envelope@plt>:
    a050:	f3 0f 1e fa          	endbr64 
    a054:	f2 ff 25 ed 29 01 00 	bnd jmp *0x129ed(%rip)        # 1ca48 <MPI_Type_get_envelope@Base>
    a05b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a060 <MPI_Alltoall@plt>:
    a060:	f3 0f 1e fa          	endbr64 
    a064:	f2 ff 25 e5 29 01 00 	bnd jmp *0x129e5(%rip)        # 1ca50 <MPI_Alltoall@Base>
    a06b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a070 <MPI_Cancel@plt>:
    a070:	f3 0f 1e fa          	endbr64 
    a074:	f2 ff 25 dd 29 01 00 	bnd jmp *0x129dd(%rip)        # 1ca58 <MPI_Cancel@Base>
    a07b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a080 <MPI_Initialized@plt>:
    a080:	f3 0f 1e fa          	endbr64 
    a084:	f2 ff 25 d5 29 01 00 	bnd jmp *0x129d5(%rip)        # 1ca60 <MPI_Initialized@Base>
    a08b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a090 <MPI_Scatterv@plt>:
    a090:	f3 0f 1e fa          	endbr64 
    a094:	f2 ff 25 cd 29 01 00 	bnd jmp *0x129cd(%rip)        # 1ca68 <MPI_Scatterv@Base>
    a09b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a0a0 <MPI_Issend@plt>:
    a0a0:	f3 0f 1e fa          	endbr64 
    a0a4:	f2 ff 25 c5 29 01 00 	bnd jmp *0x129c5(%rip)        # 1ca70 <MPI_Issend@Base>
    a0ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a0b0 <_ZN3MPI4CommC2Ev@plt>:
    a0b0:	f3 0f 1e fa          	endbr64 
    a0b4:	f2 ff 25 bd 29 01 00 	bnd jmp *0x129bd(%rip)        # 1ca78 <_ZN3MPI4CommC2Ev@Base>
    a0bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a0c0 <MPI_Graphdims_get@plt>:
    a0c0:	f3 0f 1e fa          	endbr64 
    a0c4:	f2 ff 25 b5 29 01 00 	bnd jmp *0x129b5(%rip)        # 1ca80 <MPI_Graphdims_get@Base>
    a0cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a0d0 <MPI_Info_get@plt>:
    a0d0:	f3 0f 1e fa          	endbr64 
    a0d4:	f2 ff 25 ad 29 01 00 	bnd jmp *0x129ad(%rip)        # 1ca88 <MPI_Info_get@Base>
    a0db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a0e0 <MPI_Cart_shift@plt>:
    a0e0:	f3 0f 1e fa          	endbr64 
    a0e4:	f2 ff 25 a5 29 01 00 	bnd jmp *0x129a5(%rip)        # 1ca90 <MPI_Cart_shift@Base>
    a0eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a0f0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@plt>:
    a0f0:	f3 0f 1e fa          	endbr64 
    a0f4:	f2 ff 25 9d 29 01 00 	bnd jmp *0x1299d(%rip)        # 1ca98 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@GLIBCXX_3.4.21>
    a0fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a100 <MPI_Type_get_true_extent@plt>:
    a100:	f3 0f 1e fa          	endbr64 
    a104:	f2 ff 25 95 29 01 00 	bnd jmp *0x12995(%rip)        # 1caa0 <MPI_Type_get_true_extent@Base>
    a10b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a110 <MPI_Get@plt>:
    a110:	f3 0f 1e fa          	endbr64 
    a114:	f2 ff 25 8d 29 01 00 	bnd jmp *0x1298d(%rip)        # 1caa8 <MPI_Get@Base>
    a11b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>:
    a120:	f3 0f 1e fa          	endbr64 
    a124:	f2 ff 25 85 29 01 00 	bnd jmp *0x12985(%rip)        # 1cab0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@GLIBCXX_3.4.21>
    a12b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a130 <MPI_Type_contiguous@plt>:
    a130:	f3 0f 1e fa          	endbr64 
    a134:	f2 ff 25 7d 29 01 00 	bnd jmp *0x1297d(%rip)        # 1cab8 <MPI_Type_contiguous@Base>
    a13b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a140 <MPI_Cart_sub@plt>:
    a140:	f3 0f 1e fa          	endbr64 
    a144:	f2 ff 25 75 29 01 00 	bnd jmp *0x12975(%rip)        # 1cac0 <MPI_Cart_sub@Base>
    a14b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a150 <__assert_fail@plt>:
    a150:	f3 0f 1e fa          	endbr64 
    a154:	f2 ff 25 6d 29 01 00 	bnd jmp *0x1296d(%rip)        # 1cac8 <__assert_fail@GLIBC_2.2.5>
    a15b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a160 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EmcRKS3_@plt>:
    a160:	f3 0f 1e fa          	endbr64 
    a164:	f2 ff 25 65 29 01 00 	bnd jmp *0x12965(%rip)        # 1cad0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EmcRKS3_@GLIBCXX_3.4.21>
    a16b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a170 <MPI_Scatter@plt>:
    a170:	f3 0f 1e fa          	endbr64 
    a174:	f2 ff 25 5d 29 01 00 	bnd jmp *0x1295d(%rip)        # 1cad8 <MPI_Scatter@Base>
    a17b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a180 <MPI_Win_get_group@plt>:
    a180:	f3 0f 1e fa          	endbr64 
    a184:	f2 ff 25 55 29 01 00 	bnd jmp *0x12955(%rip)        # 1cae0 <MPI_Win_get_group@Base>
    a18b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a190 <_ZNSo5flushEv@plt>:
    a190:	f3 0f 1e fa          	endbr64 
    a194:	f2 ff 25 4d 29 01 00 	bnd jmp *0x1294d(%rip)        # 1cae8 <_ZNSo5flushEv@GLIBCXX_3.4>
    a19b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a1a0 <MPI_Ssend_init@plt>:
    a1a0:	f3 0f 1e fa          	endbr64 
    a1a4:	f2 ff 25 45 29 01 00 	bnd jmp *0x12945(%rip)        # 1caf0 <MPI_Ssend_init@Base>
    a1ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a1b0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE4sizeEv@plt>:
    a1b0:	f3 0f 1e fa          	endbr64 
    a1b4:	f2 ff 25 3d 29 01 00 	bnd jmp *0x1293d(%rip)        # 1caf8 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE4sizeEv@GLIBCXX_3.4.21>
    a1bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a1c0 <MPI_Comm_spawn@plt>:
    a1c0:	f3 0f 1e fa          	endbr64 
    a1c4:	f2 ff 25 35 29 01 00 	bnd jmp *0x12935(%rip)        # 1cb00 <MPI_Comm_spawn@Base>
    a1cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a1d0 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEC1Ev@plt>:
    a1d0:	f3 0f 1e fa          	endbr64 
    a1d4:	f2 ff 25 2d 29 01 00 	bnd jmp *0x1292d(%rip)        # 1cb08 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEC1Ev@GLIBCXX_3.4.26>
    a1db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a1e0 <MPI_Barrier@plt>:
    a1e0:	f3 0f 1e fa          	endbr64 
    a1e4:	f2 ff 25 25 29 01 00 	bnd jmp *0x12925(%rip)        # 1cb10 <MPI_Barrier@Base>
    a1eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a1f0 <MPI_Win_start@plt>:
    a1f0:	f3 0f 1e fa          	endbr64 
    a1f4:	f2 ff 25 1d 29 01 00 	bnd jmp *0x1291d(%rip)        # 1cb18 <MPI_Win_start@Base>
    a1fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a200 <MPI_Pack_external@plt>:
    a200:	f3 0f 1e fa          	endbr64 
    a204:	f2 ff 25 15 29 01 00 	bnd jmp *0x12915(%rip)        # 1cb20 <MPI_Pack_external@Base>
    a20b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a210 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEED1Ev@plt>:
    a210:	f3 0f 1e fa          	endbr64 
    a214:	f2 ff 25 0d 29 01 00 	bnd jmp *0x1290d(%rip)        # 1cb28 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEED1Ev@GLIBCXX_3.4.21>
    a21b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a220 <MPI_Rsend_init@plt>:
    a220:	f3 0f 1e fa          	endbr64 
    a224:	f2 ff 25 05 29 01 00 	bnd jmp *0x12905(%rip)        # 1cb30 <MPI_Rsend_init@Base>
    a22b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a230 <MPI_Group_incl@plt>:
    a230:	f3 0f 1e fa          	endbr64 
    a234:	f2 ff 25 fd 28 01 00 	bnd jmp *0x128fd(%rip)        # 1cb38 <MPI_Group_incl@Base>
    a23b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a240 <MPI_Bsend_init@plt>:
    a240:	f3 0f 1e fa          	endbr64 
    a244:	f2 ff 25 f5 28 01 00 	bnd jmp *0x128f5(%rip)        # 1cb40 <MPI_Bsend_init@Base>
    a24b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a250 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEE7is_openEv@plt>:
    a250:	f3 0f 1e fa          	endbr64 
    a254:	f2 ff 25 ed 28 01 00 	bnd jmp *0x128ed(%rip)        # 1cb48 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEE7is_openEv@GLIBCXX_3.4>
    a25b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a260 <MPI_Sendrecv@plt>:
    a260:	f3 0f 1e fa          	endbr64 
    a264:	f2 ff 25 e5 28 01 00 	bnd jmp *0x128e5(%rip)        # 1cb50 <MPI_Sendrecv@Base>
    a26b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a270 <MPI_Finalize@plt>:
    a270:	f3 0f 1e fa          	endbr64 
    a274:	f2 ff 25 dd 28 01 00 	bnd jmp *0x128dd(%rip)        # 1cb58 <MPI_Finalize@Base>
    a27b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a280 <MPI_Request_free@plt>:
    a280:	f3 0f 1e fa          	endbr64 
    a284:	f2 ff 25 d5 28 01 00 	bnd jmp *0x128d5(%rip)        # 1cb60 <MPI_Request_free@Base>
    a28b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a290 <MPI_Gather@plt>:
    a290:	f3 0f 1e fa          	endbr64 
    a294:	f2 ff 25 cd 28 01 00 	bnd jmp *0x128cd(%rip)        # 1cb68 <MPI_Gather@Base>
    a29b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a2a0 <MPI_Grequest_complete@plt>:
    a2a0:	f3 0f 1e fa          	endbr64 
    a2a4:	f2 ff 25 c5 28 01 00 	bnd jmp *0x128c5(%rip)        # 1cb70 <MPI_Grequest_complete@Base>
    a2ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a2b0 <MPI_Type_create_darray@plt>:
    a2b0:	f3 0f 1e fa          	endbr64 
    a2b4:	f2 ff 25 bd 28 01 00 	bnd jmp *0x128bd(%rip)        # 1cb78 <MPI_Type_create_darray@Base>
    a2bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a2c0 <MPI_Win_wait@plt>:
    a2c0:	f3 0f 1e fa          	endbr64 
    a2c4:	f2 ff 25 b5 28 01 00 	bnd jmp *0x128b5(%rip)        # 1cb80 <MPI_Win_wait@Base>
    a2cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a2d0 <MPI_Rsend@plt>:
    a2d0:	f3 0f 1e fa          	endbr64 
    a2d4:	f2 ff 25 ad 28 01 00 	bnd jmp *0x128ad(%rip)        # 1cb88 <MPI_Rsend@Base>
    a2db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a2e0 <MPI_Recv@plt>:
    a2e0:	f3 0f 1e fa          	endbr64 
    a2e4:	f2 ff 25 a5 28 01 00 	bnd jmp *0x128a5(%rip)        # 1cb90 <MPI_Recv@Base>
    a2eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a2f0 <MPI_Graph_get@plt>:
    a2f0:	f3 0f 1e fa          	endbr64 
    a2f4:	f2 ff 25 9d 28 01 00 	bnd jmp *0x1289d(%rip)        # 1cb98 <MPI_Graph_get@Base>
    a2fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a300 <MPI_Comm_group@plt>:
    a300:	f3 0f 1e fa          	endbr64 
    a304:	f2 ff 25 95 28 01 00 	bnd jmp *0x12895(%rip)        # 1cba0 <MPI_Comm_group@Base>
    a30b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a310 <MPI_Comm_create@plt>:
    a310:	f3 0f 1e fa          	endbr64 
    a314:	f2 ff 25 8d 28 01 00 	bnd jmp *0x1288d(%rip)        # 1cba8 <MPI_Comm_create@Base>
    a31b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a320 <MPI_Cart_get@plt>:
    a320:	f3 0f 1e fa          	endbr64 
    a324:	f2 ff 25 85 28 01 00 	bnd jmp *0x12885(%rip)        # 1cbb0 <MPI_Cart_get@Base>
    a32b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a330 <__cxa_atexit@plt>:
    a330:	f3 0f 1e fa          	endbr64 
    a334:	f2 ff 25 7d 28 01 00 	bnd jmp *0x1287d(%rip)        # 1cbb8 <__cxa_atexit@GLIBC_2.2.5>
    a33b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a340 <MPI_Type_get_attr@plt>:
    a340:	f3 0f 1e fa          	endbr64 
    a344:	f2 ff 25 75 28 01 00 	bnd jmp *0x12875(%rip)        # 1cbc0 <MPI_Type_get_attr@Base>
    a34b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a350 <MPI_Op_commutative@plt>:
    a350:	f3 0f 1e fa          	endbr64 
    a354:	f2 ff 25 6d 28 01 00 	bnd jmp *0x1286d(%rip)        # 1cbc8 <MPI_Op_commutative@Base>
    a35b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a360 <MPI_Graph_map@plt>:
    a360:	f3 0f 1e fa          	endbr64 
    a364:	f2 ff 25 65 28 01 00 	bnd jmp *0x12865(%rip)        # 1cbd0 <MPI_Graph_map@Base>
    a36b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a370 <MPI_Win_set_name@plt>:
    a370:	f3 0f 1e fa          	endbr64 
    a374:	f2 ff 25 5d 28 01 00 	bnd jmp *0x1285d(%rip)        # 1cbd8 <MPI_Win_set_name@Base>
    a37b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a380 <MPI_Test_cancelled@plt>:
    a380:	f3 0f 1e fa          	endbr64 
    a384:	f2 ff 25 55 28 01 00 	bnd jmp *0x12855(%rip)        # 1cbe0 <MPI_Test_cancelled@Base>
    a38b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a390 <MPI_Cart_map@plt>:
    a390:	f3 0f 1e fa          	endbr64 
    a394:	f2 ff 25 4d 28 01 00 	bnd jmp *0x1284d(%rip)        # 1cbe8 <MPI_Cart_map@Base>
    a39b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a3a0 <MPI_Exscan@plt>:
    a3a0:	f3 0f 1e fa          	endbr64 
    a3a4:	f2 ff 25 45 28 01 00 	bnd jmp *0x12845(%rip)        # 1cbf0 <MPI_Exscan@Base>
    a3ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a3b0 <MPI_Pack@plt>:
    a3b0:	f3 0f 1e fa          	endbr64 
    a3b4:	f2 ff 25 3d 28 01 00 	bnd jmp *0x1283d(%rip)        # 1cbf8 <MPI_Pack@Base>
    a3bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a3c0 <__isoc99_sscanf@plt>:
    a3c0:	f3 0f 1e fa          	endbr64 
    a3c4:	f2 ff 25 35 28 01 00 	bnd jmp *0x12835(%rip)        # 1cc00 <__isoc99_sscanf@GLIBC_2.7>
    a3cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a3d0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>:
    a3d0:	f3 0f 1e fa          	endbr64 
    a3d4:	f2 ff 25 2d 28 01 00 	bnd jmp *0x1282d(%rip)        # 1cc08 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@GLIBCXX_3.4>
    a3db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a3e0 <_Znwm@plt>:
    a3e0:	f3 0f 1e fa          	endbr64 
    a3e4:	f2 ff 25 25 28 01 00 	bnd jmp *0x12825(%rip)        # 1cc10 <_Znwm@GLIBCXX_3.4>
    a3eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a3f0 <MPI_Errhandler_free@plt>:
    a3f0:	f3 0f 1e fa          	endbr64 
    a3f4:	f2 ff 25 1d 28 01 00 	bnd jmp *0x1281d(%rip)        # 1cc18 <MPI_Errhandler_free@Base>
    a3fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a400 <_ZdlPvm@plt>:
    a400:	f3 0f 1e fa          	endbr64 
    a404:	f2 ff 25 15 28 01 00 	bnd jmp *0x12815(%rip)        # 1cc20 <_ZdlPvm@CXXABI_1.3.9>
    a40b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a410 <_ZNSolsEPFRSoS_E@plt>:
    a410:	f3 0f 1e fa          	endbr64 
    a414:	f2 ff 25 0d 28 01 00 	bnd jmp *0x1280d(%rip)        # 1cc28 <_ZNSolsEPFRSoS_E@GLIBCXX_3.4>
    a41b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a420 <MPI_Abort@plt>:
    a420:	f3 0f 1e fa          	endbr64 
    a424:	f2 ff 25 05 28 01 00 	bnd jmp *0x12805(%rip)        # 1cc30 <MPI_Abort@Base>
    a42b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a430 <MPI_Intercomm_create@plt>:
    a430:	f3 0f 1e fa          	endbr64 
    a434:	f2 ff 25 fd 27 01 00 	bnd jmp *0x127fd(%rip)        # 1cc38 <MPI_Intercomm_create@Base>
    a43b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a440 <MPI_Cartdim_get@plt>:
    a440:	f3 0f 1e fa          	endbr64 
    a444:	f2 ff 25 f5 27 01 00 	bnd jmp *0x127f5(%rip)        # 1cc40 <MPI_Cartdim_get@Base>
    a44b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a450 <MPI_Cart_coords@plt>:
    a450:	f3 0f 1e fa          	endbr64 
    a454:	f2 ff 25 ed 27 01 00 	bnd jmp *0x127ed(%rip)        # 1cc48 <MPI_Cart_coords@Base>
    a45b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a460 <ompi_op_set_cxx_callback@plt>:
    a460:	f3 0f 1e fa          	endbr64 
    a464:	f2 ff 25 e5 27 01 00 	bnd jmp *0x127e5(%rip)        # 1cc50 <ompi_op_set_cxx_callback@Base>
    a46b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a470 <_ZNSaIcED1Ev@plt>:
    a470:	f3 0f 1e fa          	endbr64 
    a474:	f2 ff 25 dd 27 01 00 	bnd jmp *0x127dd(%rip)        # 1cc58 <_ZNSaIcED1Ev@GLIBCXX_3.4>
    a47b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a480 <MPI_Get_elements@plt>:
    a480:	f3 0f 1e fa          	endbr64 
    a484:	f2 ff 25 d5 27 01 00 	bnd jmp *0x127d5(%rip)        # 1cc60 <MPI_Get_elements@Base>
    a48b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a490 <MPI_Comm_delete_attr@plt>:
    a490:	f3 0f 1e fa          	endbr64 
    a494:	f2 ff 25 cd 27 01 00 	bnd jmp *0x127cd(%rip)        # 1cc68 <MPI_Comm_delete_attr@Base>
    a49b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a4a0 <MPI_Send_init@plt>:
    a4a0:	f3 0f 1e fa          	endbr64 
    a4a4:	f2 ff 25 c5 27 01 00 	bnd jmp *0x127c5(%rip)        # 1cc70 <MPI_Send_init@Base>
    a4ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a4b0 <MPI_Info_get_nkeys@plt>:
    a4b0:	f3 0f 1e fa          	endbr64 
    a4b4:	f2 ff 25 bd 27 01 00 	bnd jmp *0x127bd(%rip)        # 1cc78 <MPI_Info_get_nkeys@Base>
    a4bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a4c0 <MPI_Bsend@plt>:
    a4c0:	f3 0f 1e fa          	endbr64 
    a4c4:	f2 ff 25 b5 27 01 00 	bnd jmp *0x127b5(%rip)        # 1cc80 <MPI_Bsend@Base>
    a4cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a4d0 <MPI_Scan@plt>:
    a4d0:	f3 0f 1e fa          	endbr64 
    a4d4:	f2 ff 25 ad 27 01 00 	bnd jmp *0x127ad(%rip)        # 1cc88 <MPI_Scan@Base>
    a4db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a4e0 <__stack_chk_fail@plt>:
    a4e0:	f3 0f 1e fa          	endbr64 
    a4e4:	f2 ff 25 a5 27 01 00 	bnd jmp *0x127a5(%rip)        # 1cc90 <__stack_chk_fail@GLIBC_2.4>
    a4eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a4f0 <MPI_Comm_size@plt>:
    a4f0:	f3 0f 1e fa          	endbr64 
    a4f4:	f2 ff 25 9d 27 01 00 	bnd jmp *0x1279d(%rip)        # 1cc98 <MPI_Comm_size@Base>
    a4fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a500 <free@plt>:
    a500:	f3 0f 1e fa          	endbr64 
    a504:	f2 ff 25 95 27 01 00 	bnd jmp *0x12795(%rip)        # 1cca0 <free@GLIBC_2.2.5>
    a50b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a510 <MPI_Type_create_hindexed@plt>:
    a510:	f3 0f 1e fa          	endbr64 
    a514:	f2 ff 25 8d 27 01 00 	bnd jmp *0x1278d(%rip)        # 1cca8 <MPI_Type_create_hindexed@Base>
    a51b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a520 <MPI_Reduce_local@plt>:
    a520:	f3 0f 1e fa          	endbr64 
    a524:	f2 ff 25 85 27 01 00 	bnd jmp *0x12785(%rip)        # 1ccb0 <MPI_Reduce_local@Base>
    a52b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a530 <MPI_Reduce_scatter@plt>:
    a530:	f3 0f 1e fa          	endbr64 
    a534:	f2 ff 25 7d 27 01 00 	bnd jmp *0x1277d(%rip)        # 1ccb8 <MPI_Reduce_scatter@Base>
    a53b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a540 <exit@plt>:
    a540:	f3 0f 1e fa          	endbr64 
    a544:	f2 ff 25 75 27 01 00 	bnd jmp *0x12775(%rip)        # 1ccc0 <exit@GLIBC_2.2.5>
    a54b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a550 <MPI_Win_lock@plt>:
    a550:	f3 0f 1e fa          	endbr64 
    a554:	f2 ff 25 6d 27 01 00 	bnd jmp *0x1276d(%rip)        # 1ccc8 <MPI_Win_lock@Base>
    a55b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a560 <MPI_Allreduce@plt>:
    a560:	f3 0f 1e fa          	endbr64 
    a564:	f2 ff 25 65 27 01 00 	bnd jmp *0x12765(%rip)        # 1ccd0 <MPI_Allreduce@Base>
    a56b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a570 <MPI_Allgatherv@plt>:
    a570:	f3 0f 1e fa          	endbr64 
    a574:	f2 ff 25 5d 27 01 00 	bnd jmp *0x1275d(%rip)        # 1ccd8 <MPI_Allgatherv@Base>
    a57b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a580 <atof@plt>:
    a580:	f3 0f 1e fa          	endbr64 
    a584:	f2 ff 25 55 27 01 00 	bnd jmp *0x12755(%rip)        # 1cce0 <atof@GLIBC_2.2.5>
    a58b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a590 <MPI_Group_range_excl@plt>:
    a590:	f3 0f 1e fa          	endbr64 
    a594:	f2 ff 25 4d 27 01 00 	bnd jmp *0x1274d(%rip)        # 1cce8 <MPI_Group_range_excl@Base>
    a59b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a5a0 <_ZNKSt9basic_iosIcSt11char_traitsIcEEcvbEv@plt>:
    a5a0:	f3 0f 1e fa          	endbr64 
    a5a4:	f2 ff 25 45 27 01 00 	bnd jmp *0x12745(%rip)        # 1ccf0 <_ZNKSt9basic_iosIcSt11char_traitsIcEEcvbEv@GLIBCXX_3.4.21>
    a5ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a5b0 <MPI_Type_commit@plt>:
    a5b0:	f3 0f 1e fa          	endbr64 
    a5b4:	f2 ff 25 3d 27 01 00 	bnd jmp *0x1273d(%rip)        # 1ccf8 <MPI_Type_commit@Base>
    a5bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a5c0 <MPI_Put@plt>:
    a5c0:	f3 0f 1e fa          	endbr64 
    a5c4:	f2 ff 25 35 27 01 00 	bnd jmp *0x12735(%rip)        # 1cd00 <MPI_Put@Base>
    a5cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a5d0 <_ZdaPv@plt>:
    a5d0:	f3 0f 1e fa          	endbr64 
    a5d4:	f2 ff 25 2d 27 01 00 	bnd jmp *0x1272d(%rip)        # 1cd08 <_ZdaPv@GLIBCXX_3.4>
    a5db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a5e0 <MPI_Type_size@plt>:
    a5e0:	f3 0f 1e fa          	endbr64 
    a5e4:	f2 ff 25 25 27 01 00 	bnd jmp *0x12725(%rip)        # 1cd10 <MPI_Type_size@Base>
    a5eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a5f0 <MPI_Status_set_elements@plt>:
    a5f0:	f3 0f 1e fa          	endbr64 
    a5f4:	f2 ff 25 1d 27 01 00 	bnd jmp *0x1271d(%rip)        # 1cd18 <MPI_Status_set_elements@Base>
    a5fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a600 <MPI_Comm_get_errhandler@plt>:
    a600:	f3 0f 1e fa          	endbr64 
    a604:	f2 ff 25 15 27 01 00 	bnd jmp *0x12715(%rip)        # 1cd20 <MPI_Comm_get_errhandler@Base>
    a60b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a610 <MPI_Reduce@plt>:
    a610:	f3 0f 1e fa          	endbr64 
    a614:	f2 ff 25 0d 27 01 00 	bnd jmp *0x1270d(%rip)        # 1cd28 <MPI_Reduce@Base>
    a61b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a620 <MPI_Comm_remote_group@plt>:
    a620:	f3 0f 1e fa          	endbr64 
    a624:	f2 ff 25 05 27 01 00 	bnd jmp *0x12705(%rip)        # 1cd30 <MPI_Comm_remote_group@Base>
    a62b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a630 <malloc@plt>:
    a630:	f3 0f 1e fa          	endbr64 
    a634:	f2 ff 25 fd 26 01 00 	bnd jmp *0x126fd(%rip)        # 1cd38 <malloc@GLIBC_2.2.5>
    a63b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a640 <MPI_Wait@plt>:
    a640:	f3 0f 1e fa          	endbr64 
    a644:	f2 ff 25 f5 26 01 00 	bnd jmp *0x126f5(%rip)        # 1cd40 <MPI_Wait@Base>
    a64b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a650 <MPI_Comm_dup@plt>:
    a650:	f3 0f 1e fa          	endbr64 
    a654:	f2 ff 25 ed 26 01 00 	bnd jmp *0x126ed(%rip)        # 1cd48 <MPI_Comm_dup@Base>
    a65b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a660 <MPI_Sendrecv_replace@plt>:
    a660:	f3 0f 1e fa          	endbr64 
    a664:	f2 ff 25 e5 26 01 00 	bnd jmp *0x126e5(%rip)        # 1cd50 <MPI_Sendrecv_replace@Base>
    a66b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a670 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EPKcRKS3_@plt>:
    a670:	f3 0f 1e fa          	endbr64 
    a674:	f2 ff 25 dd 26 01 00 	bnd jmp *0x126dd(%rip)        # 1cd58 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EPKcRKS3_@GLIBCXX_3.4.21>
    a67b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a680 <MPI_Irecv@plt>:
    a680:	f3 0f 1e fa          	endbr64 
    a684:	f2 ff 25 d5 26 01 00 	bnd jmp *0x126d5(%rip)        # 1cd60 <MPI_Irecv@Base>
    a68b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a690 <MPI_Alltoallw@plt>:
    a690:	f3 0f 1e fa          	endbr64 
    a694:	f2 ff 25 cd 26 01 00 	bnd jmp *0x126cd(%rip)        # 1cd68 <MPI_Alltoallw@Base>
    a69b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a6a0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1Ev@plt>:
    a6a0:	f3 0f 1e fa          	endbr64 
    a6a4:	f2 ff 25 c5 26 01 00 	bnd jmp *0x126c5(%rip)        # 1cd70 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1Ev@GLIBCXX_3.4.21>
    a6ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a6b0 <MPI_Type_get_name@plt>:
    a6b0:	f3 0f 1e fa          	endbr64 
    a6b4:	f2 ff 25 bd 26 01 00 	bnd jmp *0x126bd(%rip)        # 1cd78 <MPI_Type_get_name@Base>
    a6bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a6c0 <MPI_Type_create_indexed_block@plt>:
    a6c0:	f3 0f 1e fa          	endbr64 
    a6c4:	f2 ff 25 b5 26 01 00 	bnd jmp *0x126b5(%rip)        # 1cd80 <MPI_Type_create_indexed_block@Base>
    a6cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a6d0 <MPI_Type_get_contents@plt>:
    a6d0:	f3 0f 1e fa          	endbr64 
    a6d4:	f2 ff 25 ad 26 01 00 	bnd jmp *0x126ad(%rip)        # 1cd88 <MPI_Type_get_contents@Base>
    a6db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a6e0 <MPI_Comm_get_attr@plt>:
    a6e0:	f3 0f 1e fa          	endbr64 
    a6e4:	f2 ff 25 a5 26 01 00 	bnd jmp *0x126a5(%rip)        # 1cd90 <MPI_Comm_get_attr@Base>
    a6eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a6f0 <__cxa_throw_bad_array_new_length@plt>:
    a6f0:	f3 0f 1e fa          	endbr64 
    a6f4:	f2 ff 25 9d 26 01 00 	bnd jmp *0x1269d(%rip)        # 1cd98 <__cxa_throw_bad_array_new_length@CXXABI_1.3.8>
    a6fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a700 <MPI_Comm_rank@plt>:
    a700:	f3 0f 1e fa          	endbr64 
    a704:	f2 ff 25 95 26 01 00 	bnd jmp *0x12695(%rip)        # 1cda0 <MPI_Comm_rank@Base>
    a70b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a710 <MPI_Init@plt>:
    a710:	f3 0f 1e fa          	endbr64 
    a714:	f2 ff 25 8d 26 01 00 	bnd jmp *0x1268d(%rip)        # 1cda8 <MPI_Init@Base>
    a71b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a720 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEE5closeEv@plt>:
    a720:	f3 0f 1e fa          	endbr64 
    a724:	f2 ff 25 85 26 01 00 	bnd jmp *0x12685(%rip)        # 1cdb0 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEE5closeEv@GLIBCXX_3.4>
    a72b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a730 <MPI_Info_delete@plt>:
    a730:	f3 0f 1e fa          	endbr64 
    a734:	f2 ff 25 7d 26 01 00 	bnd jmp *0x1267d(%rip)        # 1cdb8 <MPI_Info_delete@Base>
    a73b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a740 <MPI_Topo_test@plt>:
    a740:	f3 0f 1e fa          	endbr64 
    a744:	f2 ff 25 75 26 01 00 	bnd jmp *0x12675(%rip)        # 1cdc0 <MPI_Topo_test@Base>
    a74b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a750 <MPI_Get_count@plt>:
    a750:	f3 0f 1e fa          	endbr64 
    a754:	f2 ff 25 6d 26 01 00 	bnd jmp *0x1266d(%rip)        # 1cdc8 <MPI_Get_count@Base>
    a75b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a760 <fprintf@plt>:
    a760:	f3 0f 1e fa          	endbr64 
    a764:	f2 ff 25 65 26 01 00 	bnd jmp *0x12665(%rip)        # 1cdd0 <fprintf@GLIBC_2.2.5>
    a76b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a770 <MPI_Info_get_nthkey@plt>:
    a770:	f3 0f 1e fa          	endbr64 
    a774:	f2 ff 25 5d 26 01 00 	bnd jmp *0x1265d(%rip)        # 1cdd8 <MPI_Info_get_nthkey@Base>
    a77b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a780 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEED1Ev@plt>:
    a780:	f3 0f 1e fa          	endbr64 
    a784:	f2 ff 25 55 26 01 00 	bnd jmp *0x12655(%rip)        # 1cde0 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEED1Ev@GLIBCXX_3.4>
    a78b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a790 <MPI_Win_set_errhandler@plt>:
    a790:	f3 0f 1e fa          	endbr64 
    a794:	f2 ff 25 4d 26 01 00 	bnd jmp *0x1264d(%rip)        # 1cde8 <MPI_Win_set_errhandler@Base>
    a79b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a7a0 <MPI_Comm_remote_size@plt>:
    a7a0:	f3 0f 1e fa          	endbr64 
    a7a4:	f2 ff 25 45 26 01 00 	bnd jmp *0x12645(%rip)        # 1cdf0 <MPI_Comm_remote_size@Base>
    a7ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a7b0 <GOMP_parallel@plt>:
    a7b0:	f3 0f 1e fa          	endbr64 
    a7b4:	f2 ff 25 3d 26 01 00 	bnd jmp *0x1263d(%rip)        # 1cdf8 <GOMP_parallel@GOMP_4.0>
    a7bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a7c0 <MPI_Type_dup@plt>:
    a7c0:	f3 0f 1e fa          	endbr64 
    a7c4:	f2 ff 25 35 26 01 00 	bnd jmp *0x12635(%rip)        # 1ce00 <MPI_Type_dup@Base>
    a7cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a7d0 <_ZNSt8ios_base4InitC1Ev@plt>:
    a7d0:	f3 0f 1e fa          	endbr64 
    a7d4:	f2 ff 25 2d 26 01 00 	bnd jmp *0x1262d(%rip)        # 1ce08 <_ZNSt8ios_base4InitC1Ev@GLIBCXX_3.4>
    a7db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a7e0 <puts@plt>:
    a7e0:	f3 0f 1e fa          	endbr64 
    a7e4:	f2 ff 25 25 26 01 00 	bnd jmp *0x12625(%rip)        # 1ce10 <puts@GLIBC_2.2.5>
    a7eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a7f0 <MPI_Comm_set_errhandler@plt>:
    a7f0:	f3 0f 1e fa          	endbr64 
    a7f4:	f2 ff 25 1d 26 01 00 	bnd jmp *0x1261d(%rip)        # 1ce18 <MPI_Comm_set_errhandler@Base>
    a7fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a800 <MPI_Ibsend@plt>:
    a800:	f3 0f 1e fa          	endbr64 
    a804:	f2 ff 25 15 26 01 00 	bnd jmp *0x12615(%rip)        # 1ce20 <MPI_Ibsend@Base>
    a80b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a810 <MPI_Alltoallv@plt>:
    a810:	f3 0f 1e fa          	endbr64 
    a814:	f2 ff 25 0d 26 01 00 	bnd jmp *0x1260d(%rip)        # 1ce28 <MPI_Alltoallv@Base>
    a81b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a820 <MPI_Op_free@plt>:
    a820:	f3 0f 1e fa          	endbr64 
    a824:	f2 ff 25 05 26 01 00 	bnd jmp *0x12605(%rip)        # 1ce30 <MPI_Op_free@Base>
    a82b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a830 <strtok@plt>:
    a830:	f3 0f 1e fa          	endbr64 
    a834:	f2 ff 25 fd 25 01 00 	bnd jmp *0x125fd(%rip)        # 1ce38 <strtok@GLIBC_2.2.5>
    a83b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a840 <MPI_Comm_connect@plt>:
    a840:	f3 0f 1e fa          	endbr64 
    a844:	f2 ff 25 f5 25 01 00 	bnd jmp *0x125f5(%rip)        # 1ce40 <MPI_Comm_connect@Base>
    a84b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a850 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmRKS4_@plt>:
    a850:	f3 0f 1e fa          	endbr64 
    a854:	f2 ff 25 ed 25 01 00 	bnd jmp *0x125ed(%rip)        # 1ce48 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmRKS4_@GLIBCXX_3.4.21>
    a85b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a860 <MPI_Comm_free@plt>:
    a860:	f3 0f 1e fa          	endbr64 
    a864:	f2 ff 25 e5 25 01 00 	bnd jmp *0x125e5(%rip)        # 1ce50 <MPI_Comm_free@Base>
    a86b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a870 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEE7is_openEv@plt>:
    a870:	f3 0f 1e fa          	endbr64 
    a874:	f2 ff 25 dd 25 01 00 	bnd jmp *0x125dd(%rip)        # 1ce58 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEE7is_openEv@GLIBCXX_3.4>
    a87b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a880 <omp_get_thread_num@plt>:
    a880:	f3 0f 1e fa          	endbr64 
    a884:	f2 ff 25 d5 25 01 00 	bnd jmp *0x125d5(%rip)        # 1ce60 <omp_get_thread_num@OMP_1.0>
    a88b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a890 <MPI_Comm_split@plt>:
    a890:	f3 0f 1e fa          	endbr64 
    a894:	f2 ff 25 cd 25 01 00 	bnd jmp *0x125cd(%rip)        # 1ce68 <MPI_Comm_split@Base>
    a89b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a8a0 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEEC1ERKNSt7__cxx1112basic_stringIcS1_SaIcEEESt13_Ios_Openmode@plt>:
    a8a0:	f3 0f 1e fa          	endbr64 
    a8a4:	f2 ff 25 c5 25 01 00 	bnd jmp *0x125c5(%rip)        # 1ce70 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEEC1ERKNSt7__cxx1112basic_stringIcS1_SaIcEEESt13_Ios_Openmode@GLIBCXX_3.4.21>
    a8ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a8b0 <MPI_Type_create_subarray@plt>:
    a8b0:	f3 0f 1e fa          	endbr64 
    a8b4:	f2 ff 25 bd 25 01 00 	bnd jmp *0x125bd(%rip)        # 1ce78 <MPI_Type_create_subarray@Base>
    a8bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a8c0 <MPI_Comm_spawn_multiple@plt>:
    a8c0:	f3 0f 1e fa          	endbr64 
    a8c4:	f2 ff 25 b5 25 01 00 	bnd jmp *0x125b5(%rip)        # 1ce80 <MPI_Comm_spawn_multiple@Base>
    a8cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a8d0 <realloc@plt>:
    a8d0:	f3 0f 1e fa          	endbr64 
    a8d4:	f2 ff 25 ad 25 01 00 	bnd jmp *0x125ad(%rip)        # 1ce88 <realloc@GLIBC_2.2.5>
    a8db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a8e0 <MPI_Test@plt>:
    a8e0:	f3 0f 1e fa          	endbr64 
    a8e4:	f2 ff 25 a5 25 01 00 	bnd jmp *0x125a5(%rip)        # 1ce90 <MPI_Test@Base>
    a8eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a8f0 <MPI_Info_get_valuelen@plt>:
    a8f0:	f3 0f 1e fa          	endbr64 
    a8f4:	f2 ff 25 9d 25 01 00 	bnd jmp *0x1259d(%rip)        # 1ce98 <MPI_Info_get_valuelen@Base>
    a8fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a900 <MPI_Win_test@plt>:
    a900:	f3 0f 1e fa          	endbr64 
    a904:	f2 ff 25 95 25 01 00 	bnd jmp *0x12595(%rip)        # 1cea0 <MPI_Win_test@Base>
    a90b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a910 <MPI_Type_vector@plt>:
    a910:	f3 0f 1e fa          	endbr64 
    a914:	f2 ff 25 8d 25 01 00 	bnd jmp *0x1258d(%rip)        # 1cea8 <MPI_Type_vector@Base>
    a91b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a920 <MPI_Start@plt>:
    a920:	f3 0f 1e fa          	endbr64 
    a924:	f2 ff 25 85 25 01 00 	bnd jmp *0x12585(%rip)        # 1ceb0 <MPI_Start@Base>
    a92b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a930 <MPI_Win_call_errhandler@plt>:
    a930:	f3 0f 1e fa          	endbr64 
    a934:	f2 ff 25 7d 25 01 00 	bnd jmp *0x1257d(%rip)        # 1ceb8 <MPI_Win_call_errhandler@Base>
    a93b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a940 <MPI_Type_indexed@plt>:
    a940:	f3 0f 1e fa          	endbr64 
    a944:	f2 ff 25 75 25 01 00 	bnd jmp *0x12575(%rip)        # 1cec0 <MPI_Type_indexed@Base>
    a94b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a950 <MPI_Comm_set_name@plt>:
    a950:	f3 0f 1e fa          	endbr64 
    a954:	f2 ff 25 6d 25 01 00 	bnd jmp *0x1256d(%rip)        # 1cec8 <MPI_Comm_set_name@Base>
    a95b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a960 <MPI_Win_fence@plt>:
    a960:	f3 0f 1e fa          	endbr64 
    a964:	f2 ff 25 65 25 01 00 	bnd jmp *0x12565(%rip)        # 1ced0 <MPI_Win_fence@Base>
    a96b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a970 <MPI_Win_complete@plt>:
    a970:	f3 0f 1e fa          	endbr64 
    a974:	f2 ff 25 5d 25 01 00 	bnd jmp *0x1255d(%rip)        # 1ced8 <MPI_Win_complete@Base>
    a97b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a980 <_ZNSolsEi@plt>:
    a980:	f3 0f 1e fa          	endbr64 
    a984:	f2 ff 25 55 25 01 00 	bnd jmp *0x12555(%rip)        # 1cee0 <_ZNSolsEi@GLIBCXX_3.4>
    a98b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a990 <MPI_Group_free@plt>:
    a990:	f3 0f 1e fa          	endbr64 
    a994:	f2 ff 25 4d 25 01 00 	bnd jmp *0x1254d(%rip)        # 1cee8 <MPI_Group_free@Base>
    a99b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a9a0 <MPI_Type_create_resized@plt>:
    a9a0:	f3 0f 1e fa          	endbr64 
    a9a4:	f2 ff 25 45 25 01 00 	bnd jmp *0x12545(%rip)        # 1cef0 <MPI_Type_create_resized@Base>
    a9ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a9b0 <_Unwind_Resume@plt>:
    a9b0:	f3 0f 1e fa          	endbr64 
    a9b4:	f2 ff 25 3d 25 01 00 	bnd jmp *0x1253d(%rip)        # 1cef8 <_Unwind_Resume@GCC_3.0>
    a9bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a9c0 <_ZNSaIcEC1Ev@plt>:
    a9c0:	f3 0f 1e fa          	endbr64 
    a9c4:	f2 ff 25 35 25 01 00 	bnd jmp *0x12535(%rip)        # 1cf00 <_ZNSaIcEC1Ev@GLIBCXX_3.4>
    a9cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a9d0 <MPI_Cart_rank@plt>:
    a9d0:	f3 0f 1e fa          	endbr64 
    a9d4:	f2 ff 25 2d 25 01 00 	bnd jmp *0x1252d(%rip)        # 1cf08 <MPI_Cart_rank@Base>
    a9db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a9e0 <MPI_Comm_accept@plt>:
    a9e0:	f3 0f 1e fa          	endbr64 
    a9e4:	f2 ff 25 25 25 01 00 	bnd jmp *0x12525(%rip)        # 1cf10 <MPI_Comm_accept@Base>
    a9eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a9f0 <MPI_Win_set_attr@plt>:
    a9f0:	f3 0f 1e fa          	endbr64 
    a9f4:	f2 ff 25 1d 25 01 00 	bnd jmp *0x1251d(%rip)        # 1cf18 <MPI_Win_set_attr@Base>
    a9fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa00 <MPI_Win_delete_attr@plt>:
    aa00:	f3 0f 1e fa          	endbr64 
    aa04:	f2 ff 25 15 25 01 00 	bnd jmp *0x12515(%rip)        # 1cf20 <MPI_Win_delete_attr@Base>
    aa0b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa10 <MPI_Type_delete_attr@plt>:
    aa10:	f3 0f 1e fa          	endbr64 
    aa14:	f2 ff 25 0d 25 01 00 	bnd jmp *0x1250d(%rip)        # 1cf28 <MPI_Type_delete_attr@Base>
    aa1b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa20 <MPI_Gatherv@plt>:
    aa20:	f3 0f 1e fa          	endbr64 
    aa24:	f2 ff 25 05 25 01 00 	bnd jmp *0x12505(%rip)        # 1cf30 <MPI_Gatherv@Base>
    aa2b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa30 <MPI_Group_rank@plt>:
    aa30:	f3 0f 1e fa          	endbr64 
    aa34:	f2 ff 25 fd 24 01 00 	bnd jmp *0x124fd(%rip)        # 1cf38 <MPI_Group_rank@Base>
    aa3b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa40 <omp_get_num_threads@plt>:
    aa40:	f3 0f 1e fa          	endbr64 
    aa44:	f2 ff 25 f5 24 01 00 	bnd jmp *0x124f5(%rip)        # 1cf40 <omp_get_num_threads@OMP_1.0>
    aa4b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa50 <MPI_Intercomm_merge@plt>:
    aa50:	f3 0f 1e fa          	endbr64 
    aa54:	f2 ff 25 ed 24 01 00 	bnd jmp *0x124ed(%rip)        # 1cf48 <MPI_Intercomm_merge@Base>
    aa5b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa60 <MPI_Win_get_errhandler@plt>:
    aa60:	f3 0f 1e fa          	endbr64 
    aa64:	f2 ff 25 e5 24 01 00 	bnd jmp *0x124e5(%rip)        # 1cf50 <MPI_Win_get_errhandler@Base>
    aa6b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa70 <MPI_Accumulate@plt>:
    aa70:	f3 0f 1e fa          	endbr64 
    aa74:	f2 ff 25 dd 24 01 00 	bnd jmp *0x124dd(%rip)        # 1cf58 <MPI_Accumulate@Base>
    aa7b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa80 <_ZNSolsEd@plt>:
    aa80:	f3 0f 1e fa          	endbr64 
    aa84:	f2 ff 25 d5 24 01 00 	bnd jmp *0x124d5(%rip)        # 1cf60 <_ZNSolsEd@GLIBCXX_3.4>
    aa8b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa90 <MPI_Type_set_attr@plt>:
    aa90:	f3 0f 1e fa          	endbr64 
    aa94:	f2 ff 25 cd 24 01 00 	bnd jmp *0x124cd(%rip)        # 1cf68 <MPI_Type_set_attr@Base>
    aa9b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aaa0 <MPI_Status_set_cancelled@plt>:
    aaa0:	f3 0f 1e fa          	endbr64 
    aaa4:	f2 ff 25 c5 24 01 00 	bnd jmp *0x124c5(%rip)        # 1cf70 <MPI_Status_set_cancelled@Base>
    aaab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aab0 <MPI_Graph_neighbors_count@plt>:
    aab0:	f3 0f 1e fa          	endbr64 
    aab4:	f2 ff 25 bd 24 01 00 	bnd jmp *0x124bd(%rip)        # 1cf78 <MPI_Graph_neighbors_count@Base>
    aabb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aac0 <MPI_Info_free@plt>:
    aac0:	f3 0f 1e fa          	endbr64 
    aac4:	f2 ff 25 b5 24 01 00 	bnd jmp *0x124b5(%rip)        # 1cf80 <MPI_Info_free@Base>
    aacb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aad0 <MPI_Pack_external_size@plt>:
    aad0:	f3 0f 1e fa          	endbr64 
    aad4:	f2 ff 25 ad 24 01 00 	bnd jmp *0x124ad(%rip)        # 1cf88 <MPI_Pack_external_size@Base>
    aadb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aae0 <MPI_Type_set_name@plt>:
    aae0:	f3 0f 1e fa          	endbr64 
    aae4:	f2 ff 25 a5 24 01 00 	bnd jmp *0x124a5(%rip)        # 1cf90 <MPI_Type_set_name@Base>
    aaeb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aaf0 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEEC1ERKNSt7__cxx1112basic_stringIcS1_SaIcEEESt13_Ios_Openmode@plt>:
    aaf0:	f3 0f 1e fa          	endbr64 
    aaf4:	f2 ff 25 9d 24 01 00 	bnd jmp *0x1249d(%rip)        # 1cf98 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEEC1ERKNSt7__cxx1112basic_stringIcS1_SaIcEEESt13_Ios_Openmode@GLIBCXX_3.4.21>
    aafb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000ab00 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendERKS4_@plt>:
    ab00:	f3 0f 1e fa          	endbr64 
    ab04:	f2 ff 25 95 24 01 00 	bnd jmp *0x12495(%rip)        # 1cfa0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendERKS4_@GLIBCXX_3.4.21>
    ab0b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000ab10 <MPI_Graph_create@plt>:
    ab10:	f3 0f 1e fa          	endbr64 
    ab14:	f2 ff 25 8d 24 01 00 	bnd jmp *0x1248d(%rip)        # 1cfa8 <MPI_Graph_create@Base>
    ab1b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000ab20 <MPI_Comm_test_inter@plt>:
    ab20:	f3 0f 1e fa          	endbr64 
    ab24:	f2 ff 25 85 24 01 00 	bnd jmp *0x12485(%rip)        # 1cfb0 <MPI_Comm_test_inter@Base>
    ab2b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000ab30 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEixEm@plt>:
    ab30:	f3 0f 1e fa          	endbr64 
    ab34:	f2 ff 25 7d 24 01 00 	bnd jmp *0x1247d(%rip)        # 1cfb8 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEixEm@GLIBCXX_3.4.21>
    ab3b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

Disassembly of section .text:

000000000000ab40 <_start>:
    ab40:	f3 0f 1e fa          	endbr64 
    ab44:	31 ed                	xor    %ebp,%ebp
    ab46:	49 89 d1             	mov    %rdx,%r9
    ab49:	5e                   	pop    %rsi
    ab4a:	48 89 e2             	mov    %rsp,%rdx
    ab4d:	48 83 e4 f0          	and    $0xfffffffffffffff0,%rsp
    ab51:	50                   	push   %rax
    ab52:	54                   	push   %rsp
    ab53:	45 31 c0             	xor    %r8d,%r8d
    ab56:	31 c9                	xor    %ecx,%ecx
    ab58:	48 8d 3d ca 00 00 00 	lea    0xca(%rip),%rdi        # ac29 <main>
    ab5f:	ff 15 73 24 01 00    	call   *0x12473(%rip)        # 1cfd8 <__libc_start_main@GLIBC_2.34>
    ab65:	f4                   	hlt    
    ab66:	66 2e 0f 1f 84 00 00 	cs nopw 0x0(%rax,%rax,1)
    ab6d:	00 00 00 

000000000000ab70 <deregister_tm_clones>:
    ab70:	48 8d 3d a1 24 01 00 	lea    0x124a1(%rip),%rdi        # 1d018 <__TMC_END__>
    ab77:	48 8d 05 9a 24 01 00 	lea    0x1249a(%rip),%rax        # 1d018 <__TMC_END__>
    ab7e:	48 39 f8             	cmp    %rdi,%rax
    ab81:	74 15                	je     ab98 <deregister_tm_clones+0x28>
    ab83:	48 8b 05 56 24 01 00 	mov    0x12456(%rip),%rax        # 1cfe0 <_ITM_deregisterTMCloneTable@Base>
    ab8a:	48 85 c0             	test   %rax,%rax
    ab8d:	74 09                	je     ab98 <deregister_tm_clones+0x28>
    ab8f:	ff e0                	jmp    *%rax
    ab91:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)
    ab98:	c3                   	ret    
    ab99:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)

000000000000aba0 <register_tm_clones>:
    aba0:	48 8d 3d 71 24 01 00 	lea    0x12471(%rip),%rdi        # 1d018 <__TMC_END__>
    aba7:	48 8d 35 6a 24 01 00 	lea    0x1246a(%rip),%rsi        # 1d018 <__TMC_END__>
    abae:	48 29 fe             	sub    %rdi,%rsi
    abb1:	48 89 f0             	mov    %rsi,%rax
    abb4:	48 c1 ee 3f          	shr    $0x3f,%rsi
    abb8:	48 c1 f8 03          	sar    $0x3,%rax
    abbc:	48 01 c6             	add    %rax,%rsi
    abbf:	48 d1 fe             	sar    %rsi
    abc2:	74 14                	je     abd8 <register_tm_clones+0x38>
    abc4:	48 8b 05 25 24 01 00 	mov    0x12425(%rip),%rax        # 1cff0 <_ITM_registerTMCloneTable@Base>
    abcb:	48 85 c0             	test   %rax,%rax
    abce:	74 08                	je     abd8 <register_tm_clones+0x38>
    abd0:	ff e0                	jmp    *%rax
    abd2:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)
    abd8:	c3                   	ret    
    abd9:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)

000000000000abe0 <__do_global_dtors_aux>:
    abe0:	f3 0f 1e fa          	endbr64 
    abe4:	80 3d 65 37 01 00 00 	cmpb   $0x0,0x13765(%rip)        # 1e350 <completed.0>
    abeb:	75 2b                	jne    ac18 <__do_global_dtors_aux+0x38>
    abed:	55                   	push   %rbp
    abee:	48 83 3d d2 23 01 00 	cmpq   $0x0,0x123d2(%rip)        # 1cfc8 <__cxa_finalize@GLIBC_2.2.5>
    abf5:	00 
    abf6:	48 89 e5             	mov    %rsp,%rbp
    abf9:	74 0c                	je     ac07 <__do_global_dtors_aux+0x27>
    abfb:	48 8b 3d 06 24 01 00 	mov    0x12406(%rip),%rdi        # 1d008 <__dso_handle>
    ac02:	e8 a9 f1 ff ff       	call   9db0 <__cxa_finalize@plt>
    ac07:	e8 64 ff ff ff       	call   ab70 <deregister_tm_clones>
    ac0c:	c6 05 3d 37 01 00 01 	movb   $0x1,0x1373d(%rip)        # 1e350 <completed.0>
    ac13:	5d                   	pop    %rbp
    ac14:	c3                   	ret    
    ac15:	0f 1f 00             	nopl   (%rax)
    ac18:	c3                   	ret    
    ac19:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)

000000000000ac20 <frame_dummy>:
    ac20:	f3 0f 1e fa          	endbr64 
    ac24:	e9 77 ff ff ff       	jmp    aba0 <register_tm_clones>

000000000000ac29 <main>:
#include "dense_scheduler.h"
#include "user_def/gravity_particle.h"
#include <string>
#include "mpi.h"

int main(int argc, char** argv) {
    ac29:	f3 0f 1e fa          	endbr64 
    ac2d:	55                   	push   %rbp
    ac2e:	48 89 e5             	mov    %rsp,%rbp
    ac31:	41 55                	push   %r13
    ac33:	41 54                	push   %r12
    ac35:	53                   	push   %rbx
    ac36:	48 81 ec 78 01 00 00 	sub    $0x178,%rsp
    ac3d:	89 bd 7c fe ff ff    	mov    %edi,-0x184(%rbp)
    ac43:	48 89 b5 70 fe ff ff 	mov    %rsi,-0x190(%rbp)
    ac4a:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ac51:	00 00 
    ac53:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
    ac57:	31 c0                	xor    %eax,%eax
    const int gridX = 3, gridY=3;
    ac59:	c7 85 90 fe ff ff 03 	movl   $0x3,-0x170(%rbp)
    ac60:	00 00 00 
    ac63:	c7 85 94 fe ff ff 03 	movl   $0x3,-0x16c(%rbp)
    ac6a:	00 00 00 
    MPI_Init(&argc, &argv);
    ac6d:	48 8d 95 70 fe ff ff 	lea    -0x190(%rbp),%rdx
    ac74:	48 8d 85 7c fe ff ff 	lea    -0x184(%rbp),%rax
    ac7b:	48 89 d6             	mov    %rdx,%rsi
    ac7e:	48 89 c7             	mov    %rax,%rdi
    ac81:	e8 8a fa ff ff       	call   a710 <MPI_Init@plt>

    topology_t topology;
    setupCommunicators(&topology, gridX, gridY);
    ac86:	48 8d 85 c0 fe ff ff 	lea    -0x140(%rbp),%rax
    ac8d:	ba 03 00 00 00       	mov    $0x3,%edx
    ac92:	be 03 00 00 00       	mov    $0x3,%esi
    ac97:	48 89 c7             	mov    %rax,%rdi
    ac9a:	e8 af 7b 00 00       	call   1284e <_Z18setupCommunicatorsP10topology_tii>
    const int& rankX = topology.rankX;
    ac9f:	48 8d 85 c0 fe ff ff 	lea    -0x140(%rbp),%rax
    aca6:	48 89 85 b0 fe ff ff 	mov    %rax,-0x150(%rbp)
    const int& rankY = topology.rankY;
    acad:	48 8d 85 c0 fe ff ff 	lea    -0x140(%rbp),%rax
    acb4:	48 83 c0 04          	add    $0x4,%rax
    acb8:	48 89 85 b8 fe ff ff 	mov    %rax,-0x148(%rbp)
    const std::string inputFilename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/data100.txt";
    acbf:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    acc6:	48 89 c7             	mov    %rax,%rdi
    acc9:	e8 f2 fc ff ff       	call   a9c0 <_ZNSaIcEC1Ev@plt>
    acce:	48 8d 95 a8 fe ff ff 	lea    -0x158(%rbp),%rdx
    acd5:	48 8d 85 f0 fe ff ff 	lea    -0x110(%rbp),%rax
    acdc:	48 8d 0d 65 b5 00 00 	lea    0xb565(%rip),%rcx        # 16248 <_ZN3MPIL11LOCK_SHAREDE+0x8>
    ace3:	48 89 ce             	mov    %rcx,%rsi
    ace6:	48 89 c7             	mov    %rax,%rdi
    ace9:	e8 82 f9 ff ff       	call   a670 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EPKcRKS3_@plt>
    acee:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    acf5:	48 89 c7             	mov    %rax,%rdi
    acf8:	e8 73 f7 ff ff       	call   a470 <_ZNSaIcED1Ev@plt>
    const std::string paddingFilename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/padding100.txt";
    acfd:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    ad04:	48 89 c7             	mov    %rax,%rdi
    ad07:	e8 b4 fc ff ff       	call   a9c0 <_ZNSaIcEC1Ev@plt>
    ad0c:	48 8d 95 a8 fe ff ff 	lea    -0x158(%rbp),%rdx
    ad13:	48 8d 85 10 ff ff ff 	lea    -0xf0(%rbp),%rax
    ad1a:	48 8d 0d 7f b5 00 00 	lea    0xb57f(%rip),%rcx        # 162a0 <_ZN3MPIL11LOCK_SHAREDE+0x60>
    ad21:	48 89 ce             	mov    %rcx,%rsi
    ad24:	48 89 c7             	mov    %rax,%rdi
    ad27:	e8 44 f9 ff ff       	call   a670 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EPKcRKS3_@plt>
    ad2c:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    ad33:	48 89 c7             	mov    %rax,%rdi
    ad36:	e8 35 f7 ff ff       	call   a470 <_ZNSaIcED1Ev@plt>
    const std::string scatterFilenameX = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/scatter100.txt.x" + std::to_string(rankX) + std::to_string(rankY);
    ad3b:	48 8b 85 b8 fe ff ff 	mov    -0x148(%rbp),%rax
    ad42:	8b 10                	mov    (%rax),%edx
    ad44:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    ad48:	89 d6                	mov    %edx,%esi
    ad4a:	48 89 c7             	mov    %rax,%rdi
    ad4d:	e8 db 08 00 00       	call   b62d <_ZNSt7__cxx119to_stringEi>
    ad52:	48 8b 85 b0 fe ff ff 	mov    -0x150(%rbp),%rax
    ad59:	8b 10                	mov    (%rax),%edx
    ad5b:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    ad62:	89 d6                	mov    %edx,%esi
    ad64:	48 89 c7             	mov    %rax,%rdi
    ad67:	e8 c1 08 00 00       	call   b62d <_ZNSt7__cxx119to_stringEi>
    ad6c:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    ad70:	48 8d 95 70 ff ff ff 	lea    -0x90(%rbp),%rdx
    ad77:	48 8d 0d 7a b5 00 00 	lea    0xb57a(%rip),%rcx        # 162f8 <_ZN3MPIL11LOCK_SHAREDE+0xb8>
    ad7e:	48 89 ce             	mov    %rcx,%rsi
    ad81:	48 89 c7             	mov    %rax,%rdi
    ad84:	e8 d8 5a 00 00       	call   10861 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEPKS5_OS8_>
    ad89:	48 8d 85 30 ff ff ff 	lea    -0xd0(%rbp),%rax
    ad90:	48 8d 55 b0          	lea    -0x50(%rbp),%rdx
    ad94:	48 8d 4d 90          	lea    -0x70(%rbp),%rcx
    ad98:	48 89 ce             	mov    %rcx,%rsi
    ad9b:	48 89 c7             	mov    %rax,%rdi
    ad9e:	e8 45 59 00 00       	call   106e8 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_>
    ada3:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    ada7:	48 89 c7             	mov    %rax,%rdi
    adaa:	e8 71 f3 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    adaf:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    adb6:	48 89 c7             	mov    %rax,%rdi
    adb9:	e8 62 f3 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    adbe:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    adc2:	48 89 c7             	mov    %rax,%rdi
    adc5:	e8 56 f3 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    const std::string scatterFilenameY = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/scatter100.txt.y" + std::to_string(rankX) + std::to_string(rankY);
    adca:	48 8b 85 b8 fe ff ff 	mov    -0x148(%rbp),%rax
    add1:	8b 10                	mov    (%rax),%edx
    add3:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    add7:	89 d6                	mov    %edx,%esi
    add9:	48 89 c7             	mov    %rax,%rdi
    addc:	e8 4c 08 00 00       	call   b62d <_ZNSt7__cxx119to_stringEi>
    ade1:	48 8b 85 b0 fe ff ff 	mov    -0x150(%rbp),%rax
    ade8:	8b 10                	mov    (%rax),%edx
    adea:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    adf1:	89 d6                	mov    %edx,%esi
    adf3:	48 89 c7             	mov    %rax,%rdi
    adf6:	e8 32 08 00 00       	call   b62d <_ZNSt7__cxx119to_stringEi>
    adfb:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    adff:	48 8d 95 70 ff ff ff 	lea    -0x90(%rbp),%rdx
    ae06:	48 8d 0d 43 b5 00 00 	lea    0xb543(%rip),%rcx        # 16350 <_ZN3MPIL11LOCK_SHAREDE+0x110>
    ae0d:	48 89 ce             	mov    %rcx,%rsi
    ae10:	48 89 c7             	mov    %rax,%rdi
    ae13:	e8 49 5a 00 00       	call   10861 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEPKS5_OS8_>
    ae18:	48 8d 85 50 ff ff ff 	lea    -0xb0(%rbp),%rax
    ae1f:	48 8d 55 b0          	lea    -0x50(%rbp),%rdx
    ae23:	48 8d 4d 90          	lea    -0x70(%rbp),%rcx
    ae27:	48 89 ce             	mov    %rcx,%rsi
    ae2a:	48 89 c7             	mov    %rax,%rdi
    ae2d:	e8 b6 58 00 00       	call   106e8 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_>
    ae32:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    ae36:	48 89 c7             	mov    %rax,%rdi
    ae39:	e8 e2 f2 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    ae3e:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    ae45:	48 89 c7             	mov    %rax,%rdi
    ae48:	e8 d3 f2 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    ae4d:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    ae51:	48 89 c7             	mov    %rax,%rdi
    ae54:	e8 c7 f2 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    const std::string gatherFilename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/gather100.txt";
    ae59:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    ae60:	48 89 c7             	mov    %rax,%rdi
    ae63:	e8 58 fb ff ff       	call   a9c0 <_ZNSaIcEC1Ev@plt>
    ae68:	48 8d 95 a8 fe ff ff 	lea    -0x158(%rbp),%rdx
    ae6f:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    ae76:	48 8d 0d 2b b5 00 00 	lea    0xb52b(%rip),%rcx        # 163a8 <_ZN3MPIL11LOCK_SHAREDE+0x168>
    ae7d:	48 89 ce             	mov    %rcx,%rsi
    ae80:	48 89 c7             	mov    %rax,%rdi
    ae83:	e8 e8 f7 ff ff       	call   a670 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EPKcRKS3_@plt>
    ae88:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    ae8f:	48 89 c7             	mov    %rax,%rdi
    ae92:	e8 d9 f5 ff ff       	call   a470 <_ZNSaIcED1Ev@plt>
    const std::string shrinkFilename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/shrink100.txt";
    ae97:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    ae9e:	48 89 c7             	mov    %rax,%rdi
    aea1:	e8 1a fb ff ff       	call   a9c0 <_ZNSaIcEC1Ev@plt>
    aea6:	48 8d 95 a8 fe ff ff 	lea    -0x158(%rbp),%rdx
    aead:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    aeb1:	48 8d 0d 48 b5 00 00 	lea    0xb548(%rip),%rcx        # 16400 <_ZN3MPIL11LOCK_SHAREDE+0x1c0>
    aeb8:	48 89 ce             	mov    %rcx,%rsi
    aebb:	48 89 c7             	mov    %rax,%rdi
    aebe:	e8 ad f7 ff ff       	call   a670 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EPKcRKS3_@plt>
    aec3:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    aeca:	48 89 c7             	mov    %rax,%rdi
    aecd:	e8 9e f5 ff ff       	call   a470 <_ZNSaIcED1Ev@plt>
    chunk_particles_t* globalParticles = nullptr;
    aed2:	48 c7 85 98 fe ff ff 	movq   $0x0,-0x168(%rbp)
    aed9:	00 00 00 00 
    
    int sizeReduceXComm=-1, sizeReduceYComm=-1, sizeScatterXComm=-1, sizeScatterYComm=-1;
    aedd:	c7 85 80 fe ff ff ff 	movl   $0xffffffff,-0x180(%rbp)
    aee4:	ff ff ff 
    aee7:	c7 85 84 fe ff ff ff 	movl   $0xffffffff,-0x17c(%rbp)
    aeee:	ff ff ff 
    aef1:	c7 85 88 fe ff ff ff 	movl   $0xffffffff,-0x178(%rbp)
    aef8:	ff ff ff 
    aefb:	c7 85 8c fe ff ff ff 	movl   $0xffffffff,-0x174(%rbp)
    af02:	ff ff ff 
    if(topology.reduceXComm != MPI_COMM_NULL)
    af05:	48 8b 85 d0 fe ff ff 	mov    -0x130(%rbp),%rax
    af0c:	48 8d 15 0d 2f 01 00 	lea    0x12f0d(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    af13:	48 39 d0             	cmp    %rdx,%rax
    af16:	74 19                	je     af31 <main+0x308>
        MPI_Comm_size(topology.reduceXComm, &sizeReduceXComm);
    af18:	48 8b 85 d0 fe ff ff 	mov    -0x130(%rbp),%rax
    af1f:	48 8d 95 80 fe ff ff 	lea    -0x180(%rbp),%rdx
    af26:	48 89 d6             	mov    %rdx,%rsi
    af29:	48 89 c7             	mov    %rax,%rdi
    af2c:	e8 bf f5 ff ff       	call   a4f0 <MPI_Comm_size@plt>
    if(topology.reduceYComm != MPI_COMM_NULL)
    af31:	48 8b 85 d8 fe ff ff 	mov    -0x128(%rbp),%rax
    af38:	48 8d 15 e1 2e 01 00 	lea    0x12ee1(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    af3f:	48 39 d0             	cmp    %rdx,%rax
    af42:	74 19                	je     af5d <main+0x334>
        MPI_Comm_size(topology.reduceYComm, &sizeReduceYComm);
    af44:	48 8b 85 d8 fe ff ff 	mov    -0x128(%rbp),%rax
    af4b:	48 8d 95 84 fe ff ff 	lea    -0x17c(%rbp),%rdx
    af52:	48 89 d6             	mov    %rdx,%rsi
    af55:	48 89 c7             	mov    %rax,%rdi
    af58:	e8 93 f5 ff ff       	call   a4f0 <MPI_Comm_size@plt>
    if(topology.scatterXComm != MPI_COMM_NULL)
    af5d:	48 8b 85 e0 fe ff ff 	mov    -0x120(%rbp),%rax
    af64:	48 8d 15 b5 2e 01 00 	lea    0x12eb5(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    af6b:	48 39 d0             	cmp    %rdx,%rax
    af6e:	74 19                	je     af89 <main+0x360>
        MPI_Comm_size(topology.scatterXComm, &sizeScatterXComm);
    af70:	48 8b 85 e0 fe ff ff 	mov    -0x120(%rbp),%rax
    af77:	48 8d 95 88 fe ff ff 	lea    -0x178(%rbp),%rdx
    af7e:	48 89 d6             	mov    %rdx,%rsi
    af81:	48 89 c7             	mov    %rax,%rdi
    af84:	e8 67 f5 ff ff       	call   a4f0 <MPI_Comm_size@plt>
    if(topology.scatterYComm != MPI_COMM_NULL)
    af89:	48 8b 85 e8 fe ff ff 	mov    -0x118(%rbp),%rax
    af90:	48 8d 15 89 2e 01 00 	lea    0x12e89(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    af97:	48 39 d0             	cmp    %rdx,%rax
    af9a:	74 19                	je     afb5 <main+0x38c>
        MPI_Comm_size(topology.scatterYComm, &sizeScatterYComm);
    af9c:	48 8b 85 e8 fe ff ff 	mov    -0x118(%rbp),%rax
    afa3:	48 8d 95 8c fe ff ff 	lea    -0x174(%rbp),%rdx
    afaa:	48 89 d6             	mov    %rdx,%rsi
    afad:	48 89 c7             	mov    %rax,%rdi
    afb0:	e8 3b f5 ff ff       	call   a4f0 <MPI_Comm_size@plt>
    std::printf(
    afb5:	8b 9d 8c fe ff ff    	mov    -0x174(%rbp),%ebx
    afbb:	4c 8b 9d e8 fe ff ff 	mov    -0x118(%rbp),%r11
    afc2:	44 8b 95 88 fe ff ff 	mov    -0x178(%rbp),%r10d
    afc9:	4c 8b 8d e0 fe ff ff 	mov    -0x120(%rbp),%r9
    afd0:	44 8b 85 84 fe ff ff 	mov    -0x17c(%rbp),%r8d
    afd7:	48 8b bd d8 fe ff ff 	mov    -0x128(%rbp),%rdi
    afde:	8b b5 80 fe ff ff    	mov    -0x180(%rbp),%esi
    afe4:	4c 8b ad d0 fe ff ff 	mov    -0x130(%rbp),%r13
    afeb:	44 8b a5 c4 fe ff ff 	mov    -0x13c(%rbp),%r12d
    aff2:	8b 8d c0 fe ff ff    	mov    -0x140(%rbp),%ecx
    aff8:	8b 95 cc fe ff ff    	mov    -0x134(%rbp),%edx
    affe:	8b 85 c8 fe ff ff    	mov    -0x138(%rbp),%eax
    b004:	48 83 ec 08          	sub    $0x8,%rsp
    b008:	53                   	push   %rbx
    b009:	41 53                	push   %r11
    b00b:	41 52                	push   %r10
    b00d:	41 51                	push   %r9
    b00f:	41 50                	push   %r8
    b011:	57                   	push   %rdi
    b012:	56                   	push   %rsi
    b013:	4d 89 e9             	mov    %r13,%r9
    b016:	45 89 e0             	mov    %r12d,%r8d
    b019:	89 c6                	mov    %eax,%esi
    b01b:	48 8d 05 36 b4 00 00 	lea    0xb436(%rip),%rax        # 16458 <_ZN3MPIL11LOCK_SHAREDE+0x218>
    b022:	48 89 c7             	mov    %rax,%rdi
    b025:	b8 00 00 00 00       	mov    $0x0,%eax
    b02a:	e8 b1 ed ff ff       	call   9de0 <printf@plt>
    b02f:	48 83 c4 40          	add    $0x40,%rsp
        "(%d, %d) => (%d, %d) (%lx:%d %lx:%d %lx:%d %lx:%d)\n", 
        topology.gridX, topology.gridY, topology.rankX, topology.rankY, 
        topology.reduceXComm, sizeReduceXComm, topology.reduceYComm, sizeReduceYComm,
        topology.scatterXComm, sizeScatterXComm, topology.scatterYComm, sizeScatterYComm);
    if(globalParticles!=nullptr)
    b033:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b03a:	48 85 c0             	test   %rax,%rax
    b03d:	74 23                	je     b062 <main+0x439>
        printf("FSFS before_gather==%lx\n", globalParticles->particles);
    b03f:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b046:	48 8b 00             	mov    (%rax),%rax
    b049:	48 89 c6             	mov    %rax,%rsi
    b04c:	48 8d 05 39 b4 00 00 	lea    0xb439(%rip),%rax        # 1648c <_ZN3MPIL11LOCK_SHAREDE+0x24c>
    b053:	48 89 c7             	mov    %rax,%rdi
    b056:	b8 00 00 00 00       	mov    $0x0,%eax
    b05b:	e8 80 ed ff ff       	call   9de0 <printf@plt>
    b060:	eb 0f                	jmp    b071 <main+0x448>
    else
        printf("FSFS before_gather=null\n");
    b062:	48 8d 05 3c b4 00 00 	lea    0xb43c(%rip),%rax        # 164a5 <_ZN3MPIL11LOCK_SHAREDE+0x265>
    b069:	48 89 c7             	mov    %rax,%rdi
    b06c:	e8 6f f7 ff ff       	call   a7e0 <puts@plt>
    if(rankX==0 && rankY==0) {
    b071:	48 8b 85 b0 fe ff ff 	mov    -0x150(%rbp),%rax
    b078:	8b 00                	mov    (%rax),%eax
    b07a:	85 c0                	test   %eax,%eax
    b07c:	0f 85 95 00 00 00    	jne    b117 <main+0x4ee>
    b082:	48 8b 85 b8 fe ff ff 	mov    -0x148(%rbp),%rax
    b089:	8b 00                	mov    (%rax),%eax
    b08b:	85 c0                	test   %eax,%eax
    b08d:	0f 85 84 00 00 00    	jne    b117 <main+0x4ee>
        load_particles(inputFilename, &globalParticles);
    b093:	48 8d 95 f0 fe ff ff 	lea    -0x110(%rbp),%rdx
    b09a:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b09e:	48 89 d6             	mov    %rdx,%rsi
    b0a1:	48 89 c7             	mov    %rax,%rdi
    b0a4:	e8 37 ef ff ff       	call   9fe0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
    b0a9:	48 8d 95 98 fe ff ff 	lea    -0x168(%rbp),%rdx
    b0b0:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b0b4:	48 89 d6             	mov    %rdx,%rsi
    b0b7:	48 89 c7             	mov    %rax,%rdi
    b0ba:	e8 bd 69 00 00       	call   11a7c <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t>
    b0bf:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b0c3:	48 89 c7             	mov    %rax,%rdi
    b0c6:	e8 55 f0 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
        padding(globalParticles, gridX*gridY);
    b0cb:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b0d2:	be 09 00 00 00       	mov    $0x9,%esi
    b0d7:	48 89 c7             	mov    %rax,%rdi
    b0da:	e8 54 7d 00 00       	call   12e33 <_Z7paddingP17chunk_particles_ti>
        store_particles(paddingFilename, globalParticles);
    b0df:	48 8b 9d 98 fe ff ff 	mov    -0x168(%rbp),%rbx
    b0e6:	48 8d 95 10 ff ff ff 	lea    -0xf0(%rbp),%rdx
    b0ed:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b0f1:	48 89 d6             	mov    %rdx,%rsi
    b0f4:	48 89 c7             	mov    %rax,%rdi
    b0f7:	e8 e4 ee ff ff       	call   9fe0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
    b0fc:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b100:	48 89 de             	mov    %rbx,%rsi
    b103:	48 89 c7             	mov    %rax,%rdi
    b106:	e8 79 6a 00 00       	call   11b84 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t>
    b10b:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b10f:	48 89 c7             	mov    %rax,%rdi
    b112:	e8 09 f0 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    }
    if(globalParticles!=nullptr)
    b117:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b11e:	48 85 c0             	test   %rax,%rax
    b121:	74 23                	je     b146 <main+0x51d>
        printf("FSKS before_gather==%lx\n", globalParticles->particles);
    b123:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b12a:	48 8b 00             	mov    (%rax),%rax
    b12d:	48 89 c6             	mov    %rax,%rsi
    b130:	48 8d 05 86 b3 00 00 	lea    0xb386(%rip),%rax        # 164bd <_ZN3MPIL11LOCK_SHAREDE+0x27d>
    b137:	48 89 c7             	mov    %rax,%rdi
    b13a:	b8 00 00 00 00       	mov    $0x0,%eax
    b13f:	e8 9c ec ff ff       	call   9de0 <printf@plt>
    b144:	eb 0f                	jmp    b155 <main+0x52c>
    else
        printf("FSKS before_gather=null\n");
    b146:	48 8d 05 89 b3 00 00 	lea    0xb389(%rip),%rax        # 164d6 <_ZN3MPIL11LOCK_SHAREDE+0x296>
    b14d:	48 89 c7             	mov    %rax,%rdi
    b150:	e8 8b f6 ff ff       	call   a7e0 <puts@plt>
    chunk_particles_t *localParticlesX=nullptr, *localParticlesY=nullptr;
    b155:	48 c7 85 a0 fe ff ff 	movq   $0x0,-0x160(%rbp)
    b15c:	00 00 00 00 
    b160:	48 c7 85 a8 fe ff ff 	movq   $0x0,-0x158(%rbp)
    b167:	00 00 00 00 
    scatter(&localParticlesX, &localParticlesY, globalParticles, &topology);
    b16b:	48 8b 95 98 fe ff ff 	mov    -0x168(%rbp),%rdx
    b172:	48 8d 8d c0 fe ff ff 	lea    -0x140(%rbp),%rcx
    b179:	48 8d b5 a8 fe ff ff 	lea    -0x158(%rbp),%rsi
    b180:	48 8d 85 a0 fe ff ff 	lea    -0x160(%rbp),%rax
    b187:	48 89 c7             	mov    %rax,%rdi
    b18a:	e8 f7 8f 00 00       	call   14186 <_Z7scatterPP17chunk_particles_tS1_S0_P10topology_t>
    store_particles(scatterFilenameX, localParticlesX);
    b18f:	48 8b 9d a0 fe ff ff 	mov    -0x160(%rbp),%rbx
    b196:	48 8d 95 30 ff ff ff 	lea    -0xd0(%rbp),%rdx
    b19d:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b1a1:	48 89 d6             	mov    %rdx,%rsi
    b1a4:	48 89 c7             	mov    %rax,%rdi
    b1a7:	e8 34 ee ff ff       	call   9fe0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
    b1ac:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b1b0:	48 89 de             	mov    %rbx,%rsi
    b1b3:	48 89 c7             	mov    %rax,%rdi
    b1b6:	e8 c9 69 00 00       	call   11b84 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t>
    b1bb:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b1bf:	48 89 c7             	mov    %rax,%rdi
    b1c2:	e8 59 ef ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    store_particles(scatterFilenameY, localParticlesY);
    b1c7:	48 8b 9d a8 fe ff ff 	mov    -0x158(%rbp),%rbx
    b1ce:	48 8d 95 50 ff ff ff 	lea    -0xb0(%rbp),%rdx
    b1d5:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b1d9:	48 89 d6             	mov    %rdx,%rsi
    b1dc:	48 89 c7             	mov    %rax,%rdi
    b1df:	e8 fc ed ff ff       	call   9fe0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
    b1e4:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b1e8:	48 89 de             	mov    %rbx,%rsi
    b1eb:	48 89 c7             	mov    %rax,%rdi
    b1ee:	e8 91 69 00 00       	call   11b84 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t>
    b1f3:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b1f7:	48 89 c7             	mov    %rax,%rdi
    b1fa:	e8 21 ef ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    if(globalParticles!=nullptr)
    b1ff:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b206:	48 85 c0             	test   %rax,%rax
    b209:	74 23                	je     b22e <main+0x605>
        printf("before_gather==%lx", globalParticles->particles);
    b20b:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b212:	48 8b 00             	mov    (%rax),%rax
    b215:	48 89 c6             	mov    %rax,%rsi
    b218:	48 8d 05 cf b2 00 00 	lea    0xb2cf(%rip),%rax        # 164ee <_ZN3MPIL11LOCK_SHAREDE+0x2ae>
    b21f:	48 89 c7             	mov    %rax,%rdi
    b222:	b8 00 00 00 00       	mov    $0x0,%eax
    b227:	e8 b4 eb ff ff       	call   9de0 <printf@plt>
    b22c:	eb 14                	jmp    b242 <main+0x619>
    else
        printf("before_gather=null");
    b22e:	48 8d 05 cc b2 00 00 	lea    0xb2cc(%rip),%rax        # 16501 <_ZN3MPIL11LOCK_SHAREDE+0x2c1>
    b235:	48 89 c7             	mov    %rax,%rdi
    b238:	b8 00 00 00 00       	mov    $0x0,%eax
    b23d:	e8 9e eb ff ff       	call   9de0 <printf@plt>
    gather(localParticlesX, localParticlesY, globalParticles, &topology);
    b242:	48 8b 95 98 fe ff ff 	mov    -0x168(%rbp),%rdx
    b249:	48 8b b5 a8 fe ff ff 	mov    -0x158(%rbp),%rsi
    b250:	48 8b 85 a0 fe ff ff 	mov    -0x160(%rbp),%rax
    b257:	48 8d 8d c0 fe ff ff 	lea    -0x140(%rbp),%rcx
    b25e:	48 89 c7             	mov    %rax,%rdi
    b261:	e8 f2 97 00 00       	call   14a58 <_Z6gatherP17chunk_particles_tS0_S0_P10topology_t>
    if(globalParticles!=nullptr)
    b266:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b26d:	48 85 c0             	test   %rax,%rax
    b270:	74 23                	je     b295 <main+0x66c>
        printf(" after_gather==%lx\n", globalParticles->particles);
    b272:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b279:	48 8b 00             	mov    (%rax),%rax
    b27c:	48 89 c6             	mov    %rax,%rsi
    b27f:	48 8d 05 8e b2 00 00 	lea    0xb28e(%rip),%rax        # 16514 <_ZN3MPIL11LOCK_SHAREDE+0x2d4>
    b286:	48 89 c7             	mov    %rax,%rdi
    b289:	b8 00 00 00 00       	mov    $0x0,%eax
    b28e:	e8 4d eb ff ff       	call   9de0 <printf@plt>
    b293:	eb 0f                	jmp    b2a4 <main+0x67b>
    else
        printf(" after_gather==nullptr\n");
    b295:	48 8d 05 8c b2 00 00 	lea    0xb28c(%rip),%rax        # 16528 <_ZN3MPIL11LOCK_SHAREDE+0x2e8>
    b29c:	48 89 c7             	mov    %rax,%rdi
    b29f:	e8 3c f5 ff ff       	call   a7e0 <puts@plt>
    if(rankX==0 && rankY==0){
    b2a4:	48 8b 85 b0 fe ff ff 	mov    -0x150(%rbp),%rax
    b2ab:	8b 00                	mov    (%rax),%eax
    b2ad:	85 c0                	test   %eax,%eax
    b2af:	0f 85 89 00 00 00    	jne    b33e <main+0x715>
    b2b5:	48 8b 85 b8 fe ff ff 	mov    -0x148(%rbp),%rax
    b2bc:	8b 00                	mov    (%rax),%eax
    b2be:	85 c0                	test   %eax,%eax
    b2c0:	75 7c                	jne    b33e <main+0x715>
        store_particles(gatherFilename, globalParticles);
    b2c2:	48 8b 9d 98 fe ff ff 	mov    -0x168(%rbp),%rbx
    b2c9:	48 8d 95 70 ff ff ff 	lea    -0x90(%rbp),%rdx
    b2d0:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b2d4:	48 89 d6             	mov    %rdx,%rsi
    b2d7:	48 89 c7             	mov    %rax,%rdi
    b2da:	e8 01 ed ff ff       	call   9fe0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
    b2df:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b2e3:	48 89 de             	mov    %rbx,%rsi
    b2e6:	48 89 c7             	mov    %rax,%rdi
    b2e9:	e8 96 68 00 00       	call   11b84 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t>
    b2ee:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b2f2:	48 89 c7             	mov    %rax,%rdi
    b2f5:	e8 26 ee ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
        shrink(globalParticles);
    b2fa:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b301:	48 89 c7             	mov    %rax,%rdi
    b304:	e8 51 7f 00 00       	call   1325a <_Z6shrinkP17chunk_particles_t>
        store_particles(shrinkFilename, globalParticles);
    b309:	48 8b 9d 98 fe ff ff 	mov    -0x168(%rbp),%rbx
    b310:	48 8d 55 90          	lea    -0x70(%rbp),%rdx
    b314:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b318:	48 89 d6             	mov    %rdx,%rsi
    b31b:	48 89 c7             	mov    %rax,%rdi
    b31e:	e8 bd ec ff ff       	call   9fe0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
    b323:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b327:	48 89 de             	mov    %rbx,%rsi
    b32a:	48 89 c7             	mov    %rax,%rdi
    b32d:	e8 52 68 00 00       	call   11b84 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t>
    b332:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b336:	48 89 c7             	mov    %rax,%rdi
    b339:	e8 e2 ed ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    }

    MPI_Finalize();
    b33e:	e8 2d ef ff ff       	call   a270 <MPI_Finalize@plt>
    return 0;
    b343:	bb 00 00 00 00       	mov    $0x0,%ebx
}
    b348:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    b34c:	48 89 c7             	mov    %rax,%rdi
    b34f:	e8 cc ed ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b354:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    b35b:	48 89 c7             	mov    %rax,%rdi
    b35e:	e8 bd ed ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b363:	48 8d 85 50 ff ff ff 	lea    -0xb0(%rbp),%rax
    b36a:	48 89 c7             	mov    %rax,%rdi
    b36d:	e8 ae ed ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b372:	48 8d 85 30 ff ff ff 	lea    -0xd0(%rbp),%rax
    b379:	48 89 c7             	mov    %rax,%rdi
    b37c:	e8 9f ed ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b381:	48 8d 85 10 ff ff ff 	lea    -0xf0(%rbp),%rax
    b388:	48 89 c7             	mov    %rax,%rdi
    b38b:	e8 90 ed ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b390:	48 8d 85 f0 fe ff ff 	lea    -0x110(%rbp),%rax
    b397:	48 89 c7             	mov    %rax,%rdi
    b39a:	e8 81 ed ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b39f:	89 d8                	mov    %ebx,%eax
    b3a1:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
    b3a5:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    b3ac:	00 00 
    b3ae:	0f 84 ff 01 00 00    	je     b5b3 <main+0x98a>
    b3b4:	e9 f5 01 00 00       	jmp    b5ae <main+0x985>
    b3b9:	f3 0f 1e fa          	endbr64 
    b3bd:	48 89 c3             	mov    %rax,%rbx
    b3c0:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    b3c7:	48 89 c7             	mov    %rax,%rdi
    b3ca:	e8 a1 f0 ff ff       	call   a470 <_ZNSaIcED1Ev@plt>
    b3cf:	48 89 d8             	mov    %rbx,%rax
    b3d2:	48 89 c7             	mov    %rax,%rdi
    b3d5:	e8 d6 f5 ff ff       	call   a9b0 <_Unwind_Resume@plt>
    b3da:	f3 0f 1e fa          	endbr64 
    b3de:	48 89 c3             	mov    %rax,%rbx
    b3e1:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    b3e8:	48 89 c7             	mov    %rax,%rdi
    b3eb:	e8 80 f0 ff ff       	call   a470 <_ZNSaIcED1Ev@plt>
    b3f0:	e9 9f 01 00 00       	jmp    b594 <main+0x96b>
    b3f5:	f3 0f 1e fa          	endbr64 
    b3f9:	48 89 c3             	mov    %rax,%rbx
    b3fc:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    b400:	48 89 c7             	mov    %rax,%rdi
    b403:	e8 18 ed ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b408:	eb 07                	jmp    b411 <main+0x7e8>
    b40a:	f3 0f 1e fa          	endbr64 
    b40e:	48 89 c3             	mov    %rax,%rbx
    b411:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    b418:	48 89 c7             	mov    %rax,%rdi
    b41b:	e8 00 ed ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b420:	eb 07                	jmp    b429 <main+0x800>
    b422:	f3 0f 1e fa          	endbr64 
    b426:	48 89 c3             	mov    %rax,%rbx
    b429:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b42d:	48 89 c7             	mov    %rax,%rdi
    b430:	e8 eb ec ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b435:	e9 4b 01 00 00       	jmp    b585 <main+0x95c>
    b43a:	f3 0f 1e fa          	endbr64 
    b43e:	48 89 c3             	mov    %rax,%rbx
    b441:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    b445:	48 89 c7             	mov    %rax,%rdi
    b448:	e8 d3 ec ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b44d:	eb 07                	jmp    b456 <main+0x82d>
    b44f:	f3 0f 1e fa          	endbr64 
    b453:	48 89 c3             	mov    %rax,%rbx
    b456:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    b45d:	48 89 c7             	mov    %rax,%rdi
    b460:	e8 bb ec ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b465:	eb 07                	jmp    b46e <main+0x845>
    b467:	f3 0f 1e fa          	endbr64 
    b46b:	48 89 c3             	mov    %rax,%rbx
    b46e:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b472:	48 89 c7             	mov    %rax,%rdi
    b475:	e8 a6 ec ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b47a:	e9 ee 00 00 00       	jmp    b56d <main+0x944>
    b47f:	f3 0f 1e fa          	endbr64 
    b483:	48 89 c3             	mov    %rax,%rbx
    b486:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    b48d:	48 89 c7             	mov    %rax,%rdi
    b490:	e8 db ef ff ff       	call   a470 <_ZNSaIcED1Ev@plt>
    b495:	e9 bb 00 00 00       	jmp    b555 <main+0x92c>
    b49a:	f3 0f 1e fa          	endbr64 
    b49e:	48 89 c3             	mov    %rax,%rbx
    b4a1:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    b4a8:	48 89 c7             	mov    %rax,%rdi
    b4ab:	e8 c0 ef ff ff       	call   a470 <_ZNSaIcED1Ev@plt>
    b4b0:	e9 91 00 00 00       	jmp    b546 <main+0x91d>
    b4b5:	f3 0f 1e fa          	endbr64 
        load_particles(inputFilename, &globalParticles);
    b4b9:	48 89 c3             	mov    %rax,%rbx
    b4bc:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b4c0:	48 89 c7             	mov    %rax,%rdi
    b4c3:	e8 58 ec ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b4c8:	eb 70                	jmp    b53a <main+0x911>
    b4ca:	f3 0f 1e fa          	endbr64 
        store_particles(paddingFilename, globalParticles);
    b4ce:	48 89 c3             	mov    %rax,%rbx
    b4d1:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b4d5:	48 89 c7             	mov    %rax,%rdi
    b4d8:	e8 43 ec ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b4dd:	eb 5b                	jmp    b53a <main+0x911>
    b4df:	f3 0f 1e fa          	endbr64 
    store_particles(scatterFilenameX, localParticlesX);
    b4e3:	48 89 c3             	mov    %rax,%rbx
    b4e6:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b4ea:	48 89 c7             	mov    %rax,%rdi
    b4ed:	e8 2e ec ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b4f2:	eb 46                	jmp    b53a <main+0x911>
    b4f4:	f3 0f 1e fa          	endbr64 
    store_particles(scatterFilenameY, localParticlesY);
    b4f8:	48 89 c3             	mov    %rax,%rbx
    b4fb:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b4ff:	48 89 c7             	mov    %rax,%rdi
    b502:	e8 19 ec ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b507:	eb 31                	jmp    b53a <main+0x911>
    b509:	f3 0f 1e fa          	endbr64 
        store_particles(gatherFilename, globalParticles);
    b50d:	48 89 c3             	mov    %rax,%rbx
    b510:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b514:	48 89 c7             	mov    %rax,%rdi
    b517:	e8 04 ec ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b51c:	eb 1c                	jmp    b53a <main+0x911>
    b51e:	f3 0f 1e fa          	endbr64 
        store_particles(shrinkFilename, globalParticles);
    b522:	48 89 c3             	mov    %rax,%rbx
    b525:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b529:	48 89 c7             	mov    %rax,%rdi
    b52c:	e8 ef eb ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b531:	eb 07                	jmp    b53a <main+0x911>
    b533:	f3 0f 1e fa          	endbr64 
}
    b537:	48 89 c3             	mov    %rax,%rbx
    b53a:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    b53e:	48 89 c7             	mov    %rax,%rdi
    b541:	e8 da eb ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b546:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    b54d:	48 89 c7             	mov    %rax,%rdi
    b550:	e8 cb eb ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b555:	48 8d 85 50 ff ff ff 	lea    -0xb0(%rbp),%rax
    b55c:	48 89 c7             	mov    %rax,%rdi
    b55f:	e8 bc eb ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b564:	eb 07                	jmp    b56d <main+0x944>
    b566:	f3 0f 1e fa          	endbr64 
    b56a:	48 89 c3             	mov    %rax,%rbx
    b56d:	48 8d 85 30 ff ff ff 	lea    -0xd0(%rbp),%rax
    b574:	48 89 c7             	mov    %rax,%rdi
    b577:	e8 a4 eb ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b57c:	eb 07                	jmp    b585 <main+0x95c>
    b57e:	f3 0f 1e fa          	endbr64 
    b582:	48 89 c3             	mov    %rax,%rbx
    b585:	48 8d 85 10 ff ff ff 	lea    -0xf0(%rbp),%rax
    b58c:	48 89 c7             	mov    %rax,%rdi
    b58f:	e8 8c eb ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b594:	48 8d 85 f0 fe ff ff 	lea    -0x110(%rbp),%rax
    b59b:	48 89 c7             	mov    %rax,%rdi
    b59e:	e8 7d eb ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b5a3:	48 89 d8             	mov    %rbx,%rax
    b5a6:	48 89 c7             	mov    %rax,%rdi
    b5a9:	e8 02 f4 ff ff       	call   a9b0 <_Unwind_Resume@plt>
    b5ae:	e8 2d ef ff ff       	call   a4e0 <__stack_chk_fail@plt>
    b5b3:	48 8d 65 e8          	lea    -0x18(%rbp),%rsp
    b5b7:	5b                   	pop    %rbx
    b5b8:	41 5c                	pop    %r12
    b5ba:	41 5d                	pop    %r13
    b5bc:	5d                   	pop    %rbp
    b5bd:	c3                   	ret    

000000000000b5be <_Z41__static_initialization_and_destruction_0ii>:
    b5be:	f3 0f 1e fa          	endbr64 
    b5c2:	55                   	push   %rbp
    b5c3:	48 89 e5             	mov    %rsp,%rbp
    b5c6:	48 83 ec 10          	sub    $0x10,%rsp
    b5ca:	89 7d fc             	mov    %edi,-0x4(%rbp)
    b5cd:	89 75 f8             	mov    %esi,-0x8(%rbp)
    b5d0:	83 7d fc 01          	cmpl   $0x1,-0x4(%rbp)
    b5d4:	75 3b                	jne    b611 <_Z41__static_initialization_and_destruction_0ii+0x53>
    b5d6:	81 7d f8 ff ff 00 00 	cmpl   $0xffff,-0x8(%rbp)
    b5dd:	75 32                	jne    b611 <_Z41__static_initialization_and_destruction_0ii+0x53>
  extern wostream wclog;	/// Linked to standard error (buffered)
#endif
  ///@}

  // For construction of filebuffers for cout, cin, cerr, clog et. al.
  static ios_base::Init __ioinit;
    b5df:	48 8d 05 6b 2d 01 00 	lea    0x12d6b(%rip),%rax        # 1e351 <_ZStL8__ioinit>
    b5e6:	48 89 c7             	mov    %rax,%rdi
    b5e9:	e8 e2 f1 ff ff       	call   a7d0 <_ZNSt8ios_base4InitC1Ev@plt>
    b5ee:	48 8d 05 13 1a 01 00 	lea    0x11a13(%rip),%rax        # 1d008 <__dso_handle>
    b5f5:	48 89 c2             	mov    %rax,%rdx
    b5f8:	48 8d 05 52 2d 01 00 	lea    0x12d52(%rip),%rax        # 1e351 <_ZStL8__ioinit>
    b5ff:	48 89 c6             	mov    %rax,%rsi
    b602:	48 8b 05 ef 19 01 00 	mov    0x119ef(%rip),%rax        # 1cff8 <_ZNSt8ios_base4InitD1Ev@GLIBCXX_3.4>
    b609:	48 89 c7             	mov    %rax,%rdi
    b60c:	e8 1f ed ff ff       	call   a330 <__cxa_atexit@plt>
    b611:	90                   	nop
    b612:	c9                   	leave  
    b613:	c3                   	ret    

000000000000b614 <_GLOBAL__sub_I_main>:
    b614:	f3 0f 1e fa          	endbr64 
    b618:	55                   	push   %rbp
    b619:	48 89 e5             	mov    %rsp,%rbp
    b61c:	be ff ff 00 00       	mov    $0xffff,%esi
    b621:	bf 01 00 00 00       	mov    $0x1,%edi
    b626:	e8 93 ff ff ff       	call   b5be <_Z41__static_initialization_and_destruction_0ii>
    b62b:	5d                   	pop    %rbp
    b62c:	c3                   	ret    

000000000000b62d <_ZNSt7__cxx119to_stringEi>:

  // DR 1261. Insufficent overloads for to_string / to_wstring

  inline string
  to_string(int __val)
  {
    b62d:	f3 0f 1e fa          	endbr64 
    b631:	55                   	push   %rbp
    b632:	48 89 e5             	mov    %rsp,%rbp
    b635:	53                   	push   %rbx
    b636:	48 83 ec 38          	sub    $0x38,%rsp
    b63a:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
    b63e:	89 75 c4             	mov    %esi,-0x3c(%rbp)
    b641:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    b648:	00 00 
    b64a:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    b64e:	31 c0                	xor    %eax,%eax
    const bool __neg = __val < 0;
    b650:	8b 45 c4             	mov    -0x3c(%rbp),%eax
    b653:	c1 e8 1f             	shr    $0x1f,%eax
    b656:	88 45 df             	mov    %al,-0x21(%rbp)
    const unsigned __uval = __neg ? (unsigned)~__val + 1u : __val;
    b659:	80 7d df 00          	cmpb   $0x0,-0x21(%rbp)
    b65d:	74 07                	je     b666 <_ZNSt7__cxx119to_stringEi+0x39>
    b65f:	8b 45 c4             	mov    -0x3c(%rbp),%eax
    b662:	f7 d8                	neg    %eax
    b664:	eb 03                	jmp    b669 <_ZNSt7__cxx119to_stringEi+0x3c>
    b666:	8b 45 c4             	mov    -0x3c(%rbp),%eax
    b669:	89 45 e0             	mov    %eax,-0x20(%rbp)
    const auto __len = __detail::__to_chars_len(__uval);
    b66c:	8b 45 e0             	mov    -0x20(%rbp),%eax
    b66f:	be 0a 00 00 00       	mov    $0xa,%esi
    b674:	89 c7                	mov    %eax,%edi
    b676:	e8 c1 00 00 00       	call   b73c <_ZNSt8__detail14__to_chars_lenIjEEjT_i>
    b67b:	89 45 e4             	mov    %eax,-0x1c(%rbp)
    string __str(__neg + __len, '-');
    b67e:	48 8d 45 de          	lea    -0x22(%rbp),%rax
    b682:	48 89 c7             	mov    %rax,%rdi
    b685:	e8 36 f3 ff ff       	call   a9c0 <_ZNSaIcEC1Ev@plt>
    b68a:	0f b6 55 df          	movzbl -0x21(%rbp),%edx
    b68e:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    b691:	01 d0                	add    %edx,%eax
    b693:	89 c6                	mov    %eax,%esi
    b695:	48 8d 55 de          	lea    -0x22(%rbp),%rdx
    b699:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    b69d:	48 89 d1             	mov    %rdx,%rcx
    b6a0:	ba 2d 00 00 00       	mov    $0x2d,%edx
    b6a5:	48 89 c7             	mov    %rax,%rdi
    b6a8:	e8 b3 ea ff ff       	call   a160 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EmcRKS3_@plt>
    b6ad:	48 8d 45 de          	lea    -0x22(%rbp),%rax
    b6b1:	48 89 c7             	mov    %rax,%rdi
    b6b4:	e8 b7 ed ff ff       	call   a470 <_ZNSaIcED1Ev@plt>
    __detail::__to_chars_10_impl(&__str[__neg], __len, __uval);
    b6b9:	0f b6 55 df          	movzbl -0x21(%rbp),%edx
    b6bd:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    b6c1:	48 89 d6             	mov    %rdx,%rsi
    b6c4:	48 89 c7             	mov    %rax,%rdi
    b6c7:	e8 64 f4 ff ff       	call   ab30 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEixEm@plt>
    b6cc:	48 89 c1             	mov    %rax,%rcx
    b6cf:	8b 55 e0             	mov    -0x20(%rbp),%edx
    b6d2:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    b6d5:	89 c6                	mov    %eax,%esi
    b6d7:	48 89 cf             	mov    %rcx,%rdi
    b6da:	e8 1d 4f 00 00       	call   105fc <_ZNSt8__detail18__to_chars_10_implIjEEvPcjT_>
    return __str;
    b6df:	90                   	nop
  }
    b6e0:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    b6e4:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    b6eb:	00 00 
    b6ed:	74 43                	je     b732 <_ZNSt7__cxx119to_stringEi+0x105>
    b6ef:	eb 3c                	jmp    b72d <_ZNSt7__cxx119to_stringEi+0x100>
    b6f1:	f3 0f 1e fa          	endbr64 
    b6f5:	48 89 c3             	mov    %rax,%rbx
    b6f8:	48 8d 45 de          	lea    -0x22(%rbp),%rax
    b6fc:	48 89 c7             	mov    %rax,%rdi
    b6ff:	e8 6c ed ff ff       	call   a470 <_ZNSaIcED1Ev@plt>
    b704:	48 89 d8             	mov    %rbx,%rax
    b707:	48 89 c7             	mov    %rax,%rdi
    b70a:	e8 a1 f2 ff ff       	call   a9b0 <_Unwind_Resume@plt>
    b70f:	f3 0f 1e fa          	endbr64 
    b713:	48 89 c3             	mov    %rax,%rbx
    b716:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    b71a:	48 89 c7             	mov    %rax,%rdi
    b71d:	e8 fe e9 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b722:	48 89 d8             	mov    %rbx,%rax
    b725:	48 89 c7             	mov    %rax,%rdi
    b728:	e8 83 f2 ff ff       	call   a9b0 <_Unwind_Resume@plt>
    b72d:	e8 ae ed ff ff       	call   a4e0 <__stack_chk_fail@plt>
    b732:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    b736:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    b73a:	c9                   	leave  
    b73b:	c3                   	ret    

000000000000b73c <_ZNSt8__detail14__to_chars_lenIjEEjT_i>:
namespace __detail
{
  // Generic implementation for arbitrary bases.
  template<typename _Tp>
    _GLIBCXX14_CONSTEXPR unsigned
    __to_chars_len(_Tp __value, int __base = 10) noexcept
    b73c:	f3 0f 1e fa          	endbr64 
    b740:	55                   	push   %rbp
    b741:	48 89 e5             	mov    %rsp,%rbp
    b744:	89 7d dc             	mov    %edi,-0x24(%rbp)
    b747:	89 75 d8             	mov    %esi,-0x28(%rbp)
    {
      static_assert(is_integral<_Tp>::value, "implementation bug");
      static_assert(is_unsigned<_Tp>::value, "implementation bug");

      unsigned __n = 1;
    b74a:	c7 45 ec 01 00 00 00 	movl   $0x1,-0x14(%rbp)
      const unsigned __b2 = __base  * __base;
    b751:	8b 45 d8             	mov    -0x28(%rbp),%eax
    b754:	0f af c0             	imul   %eax,%eax
    b757:	89 45 f0             	mov    %eax,-0x10(%rbp)
      const unsigned __b3 = __b2 * __base;
    b75a:	8b 45 d8             	mov    -0x28(%rbp),%eax
    b75d:	8b 55 f0             	mov    -0x10(%rbp),%edx
    b760:	0f af c2             	imul   %edx,%eax
    b763:	89 45 f4             	mov    %eax,-0xc(%rbp)
      const unsigned long __b4 = __b3 * __base;
    b766:	8b 45 d8             	mov    -0x28(%rbp),%eax
    b769:	0f af 45 f4          	imul   -0xc(%rbp),%eax
    b76d:	89 c0                	mov    %eax,%eax
    b76f:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
      for (;;)
	{
	  if (__value < (unsigned)__base) return __n;
    b773:	8b 45 d8             	mov    -0x28(%rbp),%eax
    b776:	39 45 dc             	cmp    %eax,-0x24(%rbp)
    b779:	73 05                	jae    b780 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x44>
    b77b:	8b 45 ec             	mov    -0x14(%rbp),%eax
    b77e:	eb 46                	jmp    b7c6 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x8a>
	  if (__value < __b2) return __n + 1;
    b780:	8b 45 dc             	mov    -0x24(%rbp),%eax
    b783:	3b 45 f0             	cmp    -0x10(%rbp),%eax
    b786:	73 08                	jae    b790 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x54>
    b788:	8b 45 ec             	mov    -0x14(%rbp),%eax
    b78b:	83 c0 01             	add    $0x1,%eax
    b78e:	eb 36                	jmp    b7c6 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x8a>
	  if (__value < __b3) return __n + 2;
    b790:	8b 45 dc             	mov    -0x24(%rbp),%eax
    b793:	3b 45 f4             	cmp    -0xc(%rbp),%eax
    b796:	73 08                	jae    b7a0 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x64>
    b798:	8b 45 ec             	mov    -0x14(%rbp),%eax
    b79b:	83 c0 02             	add    $0x2,%eax
    b79e:	eb 26                	jmp    b7c6 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x8a>
	  if (__value < __b4) return __n + 3;
    b7a0:	8b 45 dc             	mov    -0x24(%rbp),%eax
    b7a3:	48 39 45 f8          	cmp    %rax,-0x8(%rbp)
    b7a7:	76 08                	jbe    b7b1 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x75>
    b7a9:	8b 45 ec             	mov    -0x14(%rbp),%eax
    b7ac:	83 c0 03             	add    $0x3,%eax
    b7af:	eb 15                	jmp    b7c6 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x8a>
	  __value /= __b4;
    b7b1:	8b 45 dc             	mov    -0x24(%rbp),%eax
    b7b4:	ba 00 00 00 00       	mov    $0x0,%edx
    b7b9:	48 f7 75 f8          	divq   -0x8(%rbp)
    b7bd:	89 45 dc             	mov    %eax,-0x24(%rbp)
	  __n += 4;
    b7c0:	83 45 ec 04          	addl   $0x4,-0x14(%rbp)
	  if (__value < (unsigned)__base) return __n;
    b7c4:	eb ad                	jmp    b773 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x37>
	}
    }
    b7c6:	5d                   	pop    %rbp
    b7c7:	c3                   	ret    

000000000000b7c8 <_ZNKSt17integral_constantIbLb1EEcvbEv>:
    struct integral_constant
    {
      static constexpr _Tp                  value = __v;
      typedef _Tp                           value_type;
      typedef integral_constant<_Tp, __v>   type;
      constexpr operator value_type() const noexcept { return value; }
    b7c8:	f3 0f 1e fa          	endbr64 
    b7cc:	55                   	push   %rbp
    b7cd:	48 89 e5             	mov    %rsp,%rbp
    b7d0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b7d4:	b8 01 00 00 00       	mov    $0x1,%eax
    b7d9:	5d                   	pop    %rbp
    b7da:	c3                   	ret    
    b7db:	90                   	nop

000000000000b7dc <_ZN3MPI8DatatypeD1Ev>:

#else

  // construction / destruction
  inline Datatype() : mpi_datatype(MPI_DATATYPE_NULL) { }
  inline virtual ~Datatype() {}
    b7dc:	f3 0f 1e fa          	endbr64 
    b7e0:	55                   	push   %rbp
    b7e1:	48 89 e5             	mov    %rsp,%rbp
    b7e4:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b7e8:	48 8d 15 81 0c 01 00 	lea    0x10c81(%rip),%rdx        # 1c470 <_ZTVN3MPI8DatatypeE+0x10>
    b7ef:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b7f3:	48 89 10             	mov    %rdx,(%rax)
    b7f6:	90                   	nop
    b7f7:	5d                   	pop    %rbp
    b7f8:	c3                   	ret    
    b7f9:	90                   	nop

000000000000b7fa <_ZN3MPI8DatatypeD0Ev>:
    b7fa:	f3 0f 1e fa          	endbr64 
    b7fe:	55                   	push   %rbp
    b7ff:	48 89 e5             	mov    %rsp,%rbp
    b802:	48 83 ec 10          	sub    $0x10,%rsp
    b806:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b80a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b80e:	48 89 c7             	mov    %rax,%rdi
    b811:	e8 c6 ff ff ff       	call   b7dc <_ZN3MPI8DatatypeD1Ev>
    b816:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b81a:	be 10 00 00 00       	mov    $0x10,%esi
    b81f:	48 89 c7             	mov    %rax,%rdi
    b822:	e8 d9 eb ff ff       	call   a400 <_ZdlPvm@plt>
    b827:	c9                   	leave  
    b828:	c3                   	ret    
    b829:	90                   	nop

000000000000b82a <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>:
  // inter-language operability
  inline Datatype(MPI_Datatype i) : mpi_datatype(i) { }
    b82a:	f3 0f 1e fa          	endbr64 
    b82e:	55                   	push   %rbp
    b82f:	48 89 e5             	mov    %rsp,%rbp
    b832:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b836:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    b83a:	48 8d 15 2f 0c 01 00 	lea    0x10c2f(%rip),%rdx        # 1c470 <_ZTVN3MPI8DatatypeE+0x10>
    b841:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b845:	48 89 10             	mov    %rdx,(%rax)
    b848:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b84c:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    b850:	48 89 50 08          	mov    %rdx,0x8(%rax)
    b854:	90                   	nop
    b855:	5d                   	pop    %rbp
    b856:	c3                   	ret    
    b857:	90                   	nop

000000000000b858 <_ZN3MPI8DatatypeaSERKP15ompi_datatype_t>:

  inline bool operator!= (const Datatype &a) const
    { return (bool) !(*this == a); }

  // inter-language operability
  inline Datatype& operator= (const MPI_Datatype &i)
    b858:	f3 0f 1e fa          	endbr64 
    b85c:	55                   	push   %rbp
    b85d:	48 89 e5             	mov    %rsp,%rbp
    b860:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b864:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    { mpi_datatype = i; return *this; }
    b868:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    b86c:	48 8b 10             	mov    (%rax),%rdx
    b86f:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b873:	48 89 50 08          	mov    %rdx,0x8(%rax)
    b877:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b87b:	5d                   	pop    %rbp
    b87c:	c3                   	ret    
    b87d:	90                   	nop

000000000000b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>:

  inline operator MPI_Datatype () const { return mpi_datatype; }
    b87e:	f3 0f 1e fa          	endbr64 
    b882:	55                   	push   %rbp
    b883:	48 89 e5             	mov    %rsp,%rbp
    b886:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b88a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b88e:	48 8b 40 08          	mov    0x8(%rax),%rax
    b892:	5d                   	pop    %rbp
    b893:	c3                   	ret    

000000000000b894 <_ZN3MPI6StatusD1Ev>:
  // copy
  Status(const Status& data) : mpi_status(data.mpi_status) { }

  Status(const MPI_Status &i) : mpi_status(i) { }

  virtual ~Status() {}
    b894:	f3 0f 1e fa          	endbr64 
    b898:	55                   	push   %rbp
    b899:	48 89 e5             	mov    %rsp,%rbp
    b89c:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b8a0:	48 8d 15 b1 fd 00 00 	lea    0xfdb1(%rip),%rdx        # 1b658 <_ZTVN3MPI6StatusE+0x10>
    b8a7:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b8ab:	48 89 10             	mov    %rdx,(%rax)
    b8ae:	90                   	nop
    b8af:	5d                   	pop    %rbp
    b8b0:	c3                   	ret    
    b8b1:	90                   	nop

000000000000b8b2 <_ZN3MPI6StatusD0Ev>:
    b8b2:	f3 0f 1e fa          	endbr64 
    b8b6:	55                   	push   %rbp
    b8b7:	48 89 e5             	mov    %rsp,%rbp
    b8ba:	48 83 ec 10          	sub    $0x10,%rsp
    b8be:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b8c2:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b8c6:	48 89 c7             	mov    %rax,%rdi
    b8c9:	e8 c6 ff ff ff       	call   b894 <_ZN3MPI6StatusD1Ev>
    b8ce:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b8d2:	be 20 00 00 00       	mov    $0x20,%esi
    b8d7:	48 89 c7             	mov    %rax,%rdi
    b8da:	e8 21 eb ff ff       	call   a400 <_ZdlPvm@plt>
    b8df:	c9                   	leave  
    b8e0:	c3                   	ret    
    b8e1:	90                   	nop

000000000000b8e2 <_ZN3MPI6StatusaSERK20ompi_status_public_t>:
    mpi_status = data.mpi_status; return *this; }

  // comparison, don't need for status

  // inter-language operability
  Status& operator= (const MPI_Status &i) {
    b8e2:	f3 0f 1e fa          	endbr64 
    b8e6:	55                   	push   %rbp
    b8e7:	48 89 e5             	mov    %rsp,%rbp
    b8ea:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b8ee:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    mpi_status = i; return *this; }
    b8f2:	48 8b 4d f8          	mov    -0x8(%rbp),%rcx
    b8f6:	48 8b 75 f0          	mov    -0x10(%rbp),%rsi
    b8fa:	48 8b 06             	mov    (%rsi),%rax
    b8fd:	48 8b 56 08          	mov    0x8(%rsi),%rdx
    b901:	48 89 41 08          	mov    %rax,0x8(%rcx)
    b905:	48 89 51 10          	mov    %rdx,0x10(%rcx)
    b909:	48 8b 46 10          	mov    0x10(%rsi),%rax
    b90d:	48 89 41 18          	mov    %rax,0x18(%rcx)
    b911:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b915:	5d                   	pop    %rbp
    b916:	c3                   	ret    
    b917:	90                   	nop

000000000000b918 <_ZN3MPI7RequestD1Ev>:

#else

  // construction / destruction
  Request() : mpi_request(MPI_REQUEST_NULL) { }
  virtual ~Request() {}
    b918:	f3 0f 1e fa          	endbr64 
    b91c:	55                   	push   %rbp
    b91d:	48 89 e5             	mov    %rsp,%rbp
    b920:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b924:	48 8d 15 35 08 01 00 	lea    0x10835(%rip),%rdx        # 1c160 <_ZTVN3MPI7RequestE+0x10>
    b92b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b92f:	48 89 10             	mov    %rdx,(%rax)
    b932:	90                   	nop
    b933:	5d                   	pop    %rbp
    b934:	c3                   	ret    
    b935:	90                   	nop

000000000000b936 <_ZN3MPI7RequestD0Ev>:
    b936:	f3 0f 1e fa          	endbr64 
    b93a:	55                   	push   %rbp
    b93b:	48 89 e5             	mov    %rsp,%rbp
    b93e:	48 83 ec 10          	sub    $0x10,%rsp
    b942:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b946:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b94a:	48 89 c7             	mov    %rax,%rdi
    b94d:	e8 c6 ff ff ff       	call   b918 <_ZN3MPI7RequestD1Ev>
    b952:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b956:	be 10 00 00 00       	mov    $0x10,%esi
    b95b:	48 89 c7             	mov    %rax,%rdi
    b95e:	e8 9d ea ff ff       	call   a400 <_ZdlPvm@plt>
    b963:	c9                   	leave  
    b964:	c3                   	ret    
    b965:	90                   	nop

000000000000b966 <_ZN3MPI7RequestC1EP14ompi_request_t>:
  Request(MPI_Request i) : mpi_request(i) { }
    b966:	f3 0f 1e fa          	endbr64 
    b96a:	55                   	push   %rbp
    b96b:	48 89 e5             	mov    %rsp,%rbp
    b96e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b972:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    b976:	48 8d 15 e3 07 01 00 	lea    0x107e3(%rip),%rdx        # 1c160 <_ZTVN3MPI7RequestE+0x10>
    b97d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b981:	48 89 10             	mov    %rdx,(%rax)
    b984:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b988:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    b98c:	48 89 50 08          	mov    %rdx,0x8(%rax)
    b990:	90                   	nop
    b991:	5d                   	pop    %rbp
    b992:	c3                   	ret    
    b993:	90                   	nop

000000000000b994 <_ZN3MPI8PrequestC1ERKP14ompi_request_t>:
    Request::operator=(r);
    pmpi_request = r.pmpi_request; return *this; }
#else
  Prequest(const Request& p) : Request(p) { }

  Prequest(const MPI_Request &i) : Request(i) { }
    b994:	f3 0f 1e fa          	endbr64 
    b998:	55                   	push   %rbp
    b999:	48 89 e5             	mov    %rsp,%rbp
    b99c:	48 83 ec 10          	sub    $0x10,%rsp
    b9a0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b9a4:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    b9a8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b9ac:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    b9b0:	48 8b 12             	mov    (%rdx),%rdx
    b9b3:	48 89 d6             	mov    %rdx,%rsi
    b9b6:	48 89 c7             	mov    %rax,%rdi
    b9b9:	e8 a8 ff ff ff       	call   b966 <_ZN3MPI7RequestC1EP14ompi_request_t>
    b9be:	48 8d 15 33 07 01 00 	lea    0x10733(%rip),%rdx        # 1c0f8 <_ZTVN3MPI8PrequestE+0x10>
    b9c5:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b9c9:	48 89 10             	mov    %rdx,(%rax)
    b9cc:	90                   	nop
    b9cd:	c9                   	leave  
    b9ce:	c3                   	ret    
    b9cf:	90                   	nop

000000000000b9d0 <_ZN3MPI8PrequestD1Ev>:

  virtual ~Prequest() { }
    b9d0:	f3 0f 1e fa          	endbr64 
    b9d4:	55                   	push   %rbp
    b9d5:	48 89 e5             	mov    %rsp,%rbp
    b9d8:	48 83 ec 10          	sub    $0x10,%rsp
    b9dc:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b9e0:	48 8d 15 11 07 01 00 	lea    0x10711(%rip),%rdx        # 1c0f8 <_ZTVN3MPI8PrequestE+0x10>
    b9e7:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b9eb:	48 89 10             	mov    %rdx,(%rax)
    b9ee:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b9f2:	48 89 c7             	mov    %rax,%rdi
    b9f5:	e8 1e ff ff ff       	call   b918 <_ZN3MPI7RequestD1Ev>
    b9fa:	90                   	nop
    b9fb:	c9                   	leave  
    b9fc:	c3                   	ret    
    b9fd:	90                   	nop

000000000000b9fe <_ZN3MPI8PrequestD0Ev>:
    b9fe:	f3 0f 1e fa          	endbr64 
    ba02:	55                   	push   %rbp
    ba03:	48 89 e5             	mov    %rsp,%rbp
    ba06:	48 83 ec 10          	sub    $0x10,%rsp
    ba0a:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    ba0e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba12:	48 89 c7             	mov    %rax,%rdi
    ba15:	e8 b6 ff ff ff       	call   b9d0 <_ZN3MPI8PrequestD1Ev>
    ba1a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba1e:	be 10 00 00 00       	mov    $0x10,%esi
    ba23:	48 89 c7             	mov    %rax,%rdi
    ba26:	e8 d5 e9 ff ff       	call   a400 <_ZdlPvm@plt>
    ba2b:	c9                   	leave  
    ba2c:	c3                   	ret    
    ba2d:	90                   	nop

000000000000ba2e <_ZN3MPI8GrequestD1Ev>:
    typedef int Cancel_function(void *, bool);

    Grequest() {}
    Grequest(const Request& req) : Request(req) {}
    Grequest(const MPI_Request &req) : Request(req) {}
    virtual ~Grequest() {}
    ba2e:	f3 0f 1e fa          	endbr64 
    ba32:	55                   	push   %rbp
    ba33:	48 89 e5             	mov    %rsp,%rbp
    ba36:	48 83 ec 10          	sub    $0x10,%rsp
    ba3a:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    ba3e:	48 8d 15 4b 06 01 00 	lea    0x1064b(%rip),%rdx        # 1c090 <_ZTVN3MPI8GrequestE+0x10>
    ba45:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba49:	48 89 10             	mov    %rdx,(%rax)
    ba4c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba50:	48 89 c7             	mov    %rax,%rdi
    ba53:	e8 c0 fe ff ff       	call   b918 <_ZN3MPI7RequestD1Ev>
    ba58:	90                   	nop
    ba59:	c9                   	leave  
    ba5a:	c3                   	ret    
    ba5b:	90                   	nop

000000000000ba5c <_ZN3MPI8GrequestD0Ev>:
    ba5c:	f3 0f 1e fa          	endbr64 
    ba60:	55                   	push   %rbp
    ba61:	48 89 e5             	mov    %rsp,%rbp
    ba64:	48 83 ec 10          	sub    $0x10,%rsp
    ba68:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    ba6c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba70:	48 89 c7             	mov    %rax,%rdi
    ba73:	e8 b6 ff ff ff       	call   ba2e <_ZN3MPI8GrequestD1Ev>
    ba78:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba7c:	be 10 00 00 00       	mov    $0x10,%esi
    ba81:	48 89 c7             	mov    %rax,%rdi
    ba84:	e8 77 e9 ff ff       	call   a400 <_ZdlPvm@plt>
    ba89:	c9                   	leave  
    ba8a:	c3                   	ret    
    ba8b:	90                   	nop

000000000000ba8c <_ZN3MPI5GroupC1EP12ompi_group_t>:
  const PMPI::Group& pmpi() { return pmpi_group; }
#else

  // construction
  inline Group() : mpi_group(MPI_GROUP_NULL) { }
  inline Group(MPI_Group i) : mpi_group(i) { }
    ba8c:	f3 0f 1e fa          	endbr64 
    ba90:	55                   	push   %rbp
    ba91:	48 89 e5             	mov    %rsp,%rbp
    ba94:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    ba98:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    ba9c:	48 8d 15 95 fc 00 00 	lea    0xfc95(%rip),%rdx        # 1b738 <_ZTVN3MPI5GroupE+0x10>
    baa3:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    baa7:	48 89 10             	mov    %rdx,(%rax)
    baaa:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    baae:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bab2:	48 89 50 08          	mov    %rdx,0x8(%rax)
    bab6:	90                   	nop
    bab7:	5d                   	pop    %rbp
    bab8:	c3                   	ret    
    bab9:	90                   	nop

000000000000baba <_ZN3MPI5GroupD1Ev>:

  // copy
  inline Group(const Group& g) : mpi_group(g.mpi_group) { }

  inline virtual ~Group() {}
    baba:	f3 0f 1e fa          	endbr64 
    babe:	55                   	push   %rbp
    babf:	48 89 e5             	mov    %rsp,%rbp
    bac2:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bac6:	48 8d 15 6b fc 00 00 	lea    0xfc6b(%rip),%rdx        # 1b738 <_ZTVN3MPI5GroupE+0x10>
    bacd:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bad1:	48 89 10             	mov    %rdx,(%rax)
    bad4:	90                   	nop
    bad5:	5d                   	pop    %rbp
    bad6:	c3                   	ret    
    bad7:	90                   	nop

000000000000bad8 <_ZN3MPI5GroupD0Ev>:
    bad8:	f3 0f 1e fa          	endbr64 
    badc:	55                   	push   %rbp
    badd:	48 89 e5             	mov    %rsp,%rbp
    bae0:	48 83 ec 10          	sub    $0x10,%rsp
    bae4:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bae8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    baec:	48 89 c7             	mov    %rax,%rdi
    baef:	e8 c6 ff ff ff       	call   baba <_ZN3MPI5GroupD1Ev>
    baf4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    baf8:	be 10 00 00 00       	mov    $0x10,%esi
    bafd:	48 89 c7             	mov    %rax,%rdi
    bb00:	e8 fb e8 ff ff       	call   a400 <_ZdlPvm@plt>
    bb05:	c9                   	leave  
    bb06:	c3                   	ret    
    bb07:	90                   	nop

000000000000bb08 <_ZNK3MPI5GroupcvP12ompi_group_tEv>:
  inline bool operator== (const Group &a) { return (bool)(mpi_group == a.mpi_group); }
  inline bool operator!= (const Group &a) { return (bool)!(*this == a); }

  // inter-language operability
  inline Group& operator= (const MPI_Group &i) { mpi_group = i; return *this; }
  inline operator MPI_Group () const { return mpi_group; }
    bb08:	f3 0f 1e fa          	endbr64 
    bb0c:	55                   	push   %rbp
    bb0d:	48 89 e5             	mov    %rsp,%rbp
    bb10:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bb14:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bb18:	48 8b 40 08          	mov    0x8(%rax),%rax
    bb1c:	5d                   	pop    %rbp
    bb1d:	c3                   	ret    

000000000000bb1e <_ZN3MPI9Comm_NullC1EP19ompi_communicator_t>:
  // construction
  inline Comm_Null() : mpi_comm(MPI_COMM_NULL) { }
  // copy
  inline Comm_Null(const Comm_Null& data) : mpi_comm(data.mpi_comm) { }
  // inter-language operability
  inline Comm_Null(MPI_Comm data) : mpi_comm(data) { }
    bb1e:	f3 0f 1e fa          	endbr64 
    bb22:	55                   	push   %rbp
    bb23:	48 89 e5             	mov    %rsp,%rbp
    bb26:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bb2a:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    bb2e:	48 8d 15 8b 06 01 00 	lea    0x1068b(%rip),%rdx        # 1c1c0 <_ZTVN3MPI9Comm_NullE+0x10>
    bb35:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bb39:	48 89 10             	mov    %rdx,(%rax)
    bb3c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bb40:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bb44:	48 89 50 08          	mov    %rdx,0x8(%rax)
    bb48:	90                   	nop
    bb49:	5d                   	pop    %rbp
    bb4a:	c3                   	ret    
    bb4b:	90                   	nop

000000000000bb4c <_ZN3MPI9Comm_NullD1Ev>:

  // destruction
  virtual inline ~Comm_Null() { }
    bb4c:	f3 0f 1e fa          	endbr64 
    bb50:	55                   	push   %rbp
    bb51:	48 89 e5             	mov    %rsp,%rbp
    bb54:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bb58:	48 8d 15 61 06 01 00 	lea    0x10661(%rip),%rdx        # 1c1c0 <_ZTVN3MPI9Comm_NullE+0x10>
    bb5f:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bb63:	48 89 10             	mov    %rdx,(%rax)
    bb66:	90                   	nop
    bb67:	5d                   	pop    %rbp
    bb68:	c3                   	ret    
    bb69:	90                   	nop

000000000000bb6a <_ZN3MPI9Comm_NullD0Ev>:
    bb6a:	f3 0f 1e fa          	endbr64 
    bb6e:	55                   	push   %rbp
    bb6f:	48 89 e5             	mov    %rsp,%rbp
    bb72:	48 83 ec 10          	sub    $0x10,%rsp
    bb76:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bb7a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bb7e:	48 89 c7             	mov    %rax,%rdi
    bb81:	e8 c6 ff ff ff       	call   bb4c <_ZN3MPI9Comm_NullD1Ev>
    bb86:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bb8a:	be 10 00 00 00       	mov    $0x10,%esi
    bb8f:	48 89 c7             	mov    %rax,%rdi
    bb92:	e8 69 e8 ff ff       	call   a400 <_ZdlPvm@plt>
    bb97:	c9                   	leave  
    bb98:	c3                   	ret    
    bb99:	90                   	nop

000000000000bb9a <_ZNK3MPI9Comm_NullcvP19ompi_communicator_tEv>:

  inline bool operator!=(const Comm_Null& data) const {
    return (bool) !(*this == data);}

  // inter-language operability (conversion operators)
  inline operator MPI_Comm() const { return mpi_comm; }
    bb9a:	f3 0f 1e fa          	endbr64 
    bb9e:	55                   	push   %rbp
    bb9f:	48 89 e5             	mov    %rsp,%rbp
    bba2:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bba6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bbaa:	48 8b 40 08          	mov    0x8(%rax),%rax
    bbae:	5d                   	pop    %rbp
    bbaf:	c3                   	ret    

000000000000bbb0 <_ZN3MPI4CommC1EP19ompi_communicator_t>:
  }

#else
  Comm(const Comm& data) : Comm_Null(data.mpi_comm) { }
  // inter-language operability
  Comm(MPI_Comm data) : Comm_Null(data) { }
    bbb0:	f3 0f 1e fa          	endbr64 
    bbb4:	55                   	push   %rbp
    bbb5:	48 89 e5             	mov    %rsp,%rbp
    bbb8:	48 83 ec 10          	sub    $0x10,%rsp
    bbbc:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bbc0:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    bbc4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bbc8:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bbcc:	48 89 d6             	mov    %rdx,%rsi
    bbcf:	48 89 c7             	mov    %rax,%rdi
    bbd2:	e8 47 ff ff ff       	call   bb1e <_ZN3MPI9Comm_NullC1EP19ompi_communicator_t>
    bbd7:	48 8d 15 c2 06 01 00 	lea    0x106c2(%rip),%rdx        # 1c2a0 <_ZTVN3MPI4CommE+0x10>
    bbde:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bbe2:	48 89 10             	mov    %rdx,(%rax)
    bbe5:	90                   	nop
    bbe6:	c9                   	leave  
    bbe7:	c3                   	ret    

000000000000bbe8 <_ZN3MPI3WinD1Ev>:
  // copy
  Win(const Win& data) : mpi_win(data.mpi_win) { }

  Win(MPI_Win i) : mpi_win(i) { }

  virtual ~Win() { }
    bbe8:	f3 0f 1e fa          	endbr64 
    bbec:	55                   	push   %rbp
    bbed:	48 89 e5             	mov    %rsp,%rbp
    bbf0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bbf4:	48 8d 15 e5 05 01 00 	lea    0x105e5(%rip),%rdx        # 1c1e0 <_ZTVN3MPI3WinE+0x10>
    bbfb:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bbff:	48 89 10             	mov    %rdx,(%rax)
    bc02:	90                   	nop
    bc03:	5d                   	pop    %rbp
    bc04:	c3                   	ret    
    bc05:	90                   	nop

000000000000bc06 <_ZN3MPI3WinD0Ev>:
    bc06:	f3 0f 1e fa          	endbr64 
    bc0a:	55                   	push   %rbp
    bc0b:	48 89 e5             	mov    %rsp,%rbp
    bc0e:	48 83 ec 10          	sub    $0x10,%rsp
    bc12:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bc16:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bc1a:	48 89 c7             	mov    %rax,%rdi
    bc1d:	e8 c6 ff ff ff       	call   bbe8 <_ZN3MPI3WinD1Ev>
    bc22:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bc26:	be 10 00 00 00       	mov    $0x10,%esi
    bc2b:	48 89 c7             	mov    %rax,%rdi
    bc2e:	e8 cd e7 ff ff       	call   a400 <_ZdlPvm@plt>
    bc33:	c9                   	leave  
    bc34:	c3                   	ret    
    bc35:	90                   	nop

000000000000bc36 <_ZN3MPI10ErrhandlerD1Ev>:
public:
  // construction / destruction
  inline Errhandler()
    : mpi_errhandler(MPI_ERRHANDLER_NULL) {}

  inline virtual ~Errhandler() { }
    bc36:	f3 0f 1e fa          	endbr64 
    bc3a:	55                   	push   %rbp
    bc3b:	48 89 e5             	mov    %rsp,%rbp
    bc3e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bc42:	48 8d 15 87 fa 00 00 	lea    0xfa87(%rip),%rdx        # 1b6d0 <_ZTVN3MPI10ErrhandlerE+0x10>
    bc49:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bc4d:	48 89 10             	mov    %rdx,(%rax)
    bc50:	90                   	nop
    bc51:	5d                   	pop    %rbp
    bc52:	c3                   	ret    
    bc53:	90                   	nop

000000000000bc54 <_ZN3MPI10ErrhandlerD0Ev>:
    bc54:	f3 0f 1e fa          	endbr64 
    bc58:	55                   	push   %rbp
    bc59:	48 89 e5             	mov    %rsp,%rbp
    bc5c:	48 83 ec 10          	sub    $0x10,%rsp
    bc60:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bc64:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bc68:	48 89 c7             	mov    %rax,%rdi
    bc6b:	e8 c6 ff ff ff       	call   bc36 <_ZN3MPI10ErrhandlerD1Ev>
    bc70:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bc74:	be 10 00 00 00       	mov    $0x10,%esi
    bc79:	48 89 c7             	mov    %rax,%rdi
    bc7c:	e8 7f e7 ff ff       	call   a400 <_ZdlPvm@plt>
    bc81:	c9                   	leave  
    bc82:	c3                   	ret    
    bc83:	90                   	nop

000000000000bc84 <_ZN3MPI10ErrhandlerC1EP17ompi_errhandler_t>:

  inline Errhandler(MPI_Errhandler i)
    bc84:	f3 0f 1e fa          	endbr64 
    bc88:	55                   	push   %rbp
    bc89:	48 89 e5             	mov    %rsp,%rbp
    bc8c:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bc90:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    : mpi_errhandler(i) {}
    bc94:	48 8d 15 35 fa 00 00 	lea    0xfa35(%rip),%rdx        # 1b6d0 <_ZTVN3MPI10ErrhandlerE+0x10>
    bc9b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bc9f:	48 89 10             	mov    %rdx,(%rax)
    bca2:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bca6:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bcaa:	48 89 50 08          	mov    %rdx,0x8(%rax)
    bcae:	90                   	nop
    bcaf:	5d                   	pop    %rbp
    bcb0:	c3                   	ret    
    bcb1:	90                   	nop

000000000000bcb2 <_ZNK3MPI10ErrhandlercvP17ompi_errhandler_tEv>:

  // inter-language operability
  inline Errhandler& operator= (const MPI_Errhandler &i) {
    mpi_errhandler = i; return *this; }

  inline operator MPI_Errhandler() const { return mpi_errhandler; }
    bcb2:	f3 0f 1e fa          	endbr64 
    bcb6:	55                   	push   %rbp
    bcb7:	48 89 e5             	mov    %rsp,%rbp
    bcba:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bcbe:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bcc2:	48 8b 40 08          	mov    0x8(%rax),%rax
    bcc6:	5d                   	pop    %rbp
    bcc7:	c3                   	ret    

000000000000bcc8 <_ZN3MPI4CommD1Ev>:
class Comm : public Comm_Null {
    bcc8:	f3 0f 1e fa          	endbr64 
    bccc:	55                   	push   %rbp
    bccd:	48 89 e5             	mov    %rsp,%rbp
    bcd0:	48 83 ec 10          	sub    $0x10,%rsp
    bcd4:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bcd8:	48 8d 15 c1 05 01 00 	lea    0x105c1(%rip),%rdx        # 1c2a0 <_ZTVN3MPI4CommE+0x10>
    bcdf:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bce3:	48 89 10             	mov    %rdx,(%rax)
    bce6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bcea:	48 89 c7             	mov    %rax,%rdi
    bced:	e8 5a fe ff ff       	call   bb4c <_ZN3MPI9Comm_NullD1Ev>
    bcf2:	90                   	nop
    bcf3:	c9                   	leave  
    bcf4:	c3                   	ret    
    bcf5:	90                   	nop

000000000000bcf6 <_ZN3MPI4CommD0Ev>:
    bcf6:	f3 0f 1e fa          	endbr64 
    bcfa:	55                   	push   %rbp
    bcfb:	48 89 e5             	mov    %rsp,%rbp
    bcfe:	48 83 ec 10          	sub    $0x10,%rsp
    bd02:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bd06:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd0a:	48 89 c7             	mov    %rax,%rdi
    bd0d:	e8 b6 ff ff ff       	call   bcc8 <_ZN3MPI4CommD1Ev>
    bd12:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd16:	be 10 00 00 00       	mov    $0x10,%esi
    bd1b:	48 89 c7             	mov    %rax,%rdi
    bd1e:	e8 dd e6 ff ff       	call   a400 <_ZdlPvm@plt>
    bd23:	c9                   	leave  
    bd24:	c3                   	ret    
    bd25:	90                   	nop

000000000000bd26 <_ZN3MPI9IntracommC1Ev>:

class Intracomm : public Comm {
public:

  // construction
  Intracomm() { }
    bd26:	f3 0f 1e fa          	endbr64 
    bd2a:	55                   	push   %rbp
    bd2b:	48 89 e5             	mov    %rsp,%rbp
    bd2e:	48 83 ec 10          	sub    $0x10,%rsp
    bd32:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bd36:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd3a:	48 89 c7             	mov    %rax,%rdi
    bd3d:	e8 6e e3 ff ff       	call   a0b0 <_ZN3MPI4CommC2Ev@plt>
    bd42:	48 8d 15 0f 01 01 00 	lea    0x1010f(%rip),%rdx        # 1be58 <_ZTVN3MPI9IntracommE+0x10>
    bd49:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd4d:	48 89 10             	mov    %rdx,(%rax)
    bd50:	90                   	nop
    bd51:	c9                   	leave  
    bd52:	c3                   	ret    
    bd53:	90                   	nop

000000000000bd54 <_ZN3MPI9IntracommD1Ev>:
class Intracomm : public Comm {
    bd54:	f3 0f 1e fa          	endbr64 
    bd58:	55                   	push   %rbp
    bd59:	48 89 e5             	mov    %rsp,%rbp
    bd5c:	48 83 ec 10          	sub    $0x10,%rsp
    bd60:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bd64:	48 8d 15 ed 00 01 00 	lea    0x100ed(%rip),%rdx        # 1be58 <_ZTVN3MPI9IntracommE+0x10>
    bd6b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd6f:	48 89 10             	mov    %rdx,(%rax)
    bd72:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd76:	48 89 c7             	mov    %rax,%rdi
    bd79:	e8 4a ff ff ff       	call   bcc8 <_ZN3MPI4CommD1Ev>
    bd7e:	90                   	nop
    bd7f:	c9                   	leave  
    bd80:	c3                   	ret    
    bd81:	90                   	nop

000000000000bd82 <_ZN3MPI9IntracommD0Ev>:
    bd82:	f3 0f 1e fa          	endbr64 
    bd86:	55                   	push   %rbp
    bd87:	48 89 e5             	mov    %rsp,%rbp
    bd8a:	48 83 ec 10          	sub    $0x10,%rsp
    bd8e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bd92:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd96:	48 89 c7             	mov    %rax,%rdi
    bd99:	e8 b6 ff ff ff       	call   bd54 <_ZN3MPI9IntracommD1Ev>
    bd9e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bda2:	be 10 00 00 00       	mov    $0x10,%esi
    bda7:	48 89 c7             	mov    %rax,%rdi
    bdaa:	e8 51 e6 ff ff       	call   a400 <_ZdlPvm@plt>
    bdaf:	c9                   	leave  
    bdb0:	c3                   	ret    
    bdb1:	90                   	nop

000000000000bdb2 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>:
  // construction
  Intercomm() : Comm(MPI_COMM_NULL) { }
  // copy
  Intercomm(const Comm_Null& data) : Comm(data) { }
  // inter-language operability
  Intercomm(MPI_Comm data) : Comm(data) { }
    bdb2:	f3 0f 1e fa          	endbr64 
    bdb6:	55                   	push   %rbp
    bdb7:	48 89 e5             	mov    %rsp,%rbp
    bdba:	48 83 ec 10          	sub    $0x10,%rsp
    bdbe:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bdc2:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    bdc6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bdca:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bdce:	48 89 d6             	mov    %rdx,%rsi
    bdd1:	48 89 c7             	mov    %rax,%rdi
    bdd4:	e8 d7 fd ff ff       	call   bbb0 <_ZN3MPI4CommC1EP19ompi_communicator_t>
    bdd9:	48 8d 15 b0 f9 00 00 	lea    0xf9b0(%rip),%rdx        # 1b790 <_ZTVN3MPI9IntercommE+0x10>
    bde0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bde4:	48 89 10             	mov    %rdx,(%rax)
    bde7:	90                   	nop
    bde8:	c9                   	leave  
    bde9:	c3                   	ret    

000000000000bdea <_ZN3MPI4InfoC1EP11ompi_info_t>:

  Info() : mpi_info(MPI_INFO_NULL) { }
  // copy
  Info(const Info& data) : mpi_info(data.mpi_info) { }

  Info(MPI_Info i) : mpi_info(i) { }
    bdea:	f3 0f 1e fa          	endbr64 
    bdee:	55                   	push   %rbp
    bdef:	48 89 e5             	mov    %rsp,%rbp
    bdf2:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bdf6:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    bdfa:	48 8d 15 f7 f7 00 00 	lea    0xf7f7(%rip),%rdx        # 1b5f8 <_ZTVN3MPI4InfoE+0x10>
    be01:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    be05:	48 89 10             	mov    %rdx,(%rax)
    be08:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    be0c:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    be10:	48 89 50 08          	mov    %rdx,0x8(%rax)
    be14:	90                   	nop
    be15:	5d                   	pop    %rbp
    be16:	c3                   	ret    
    be17:	90                   	nop

000000000000be18 <_ZN3MPI4InfoD1Ev>:

  virtual ~Info() {}
    be18:	f3 0f 1e fa          	endbr64 
    be1c:	55                   	push   %rbp
    be1d:	48 89 e5             	mov    %rsp,%rbp
    be20:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    be24:	48 8d 15 cd f7 00 00 	lea    0xf7cd(%rip),%rdx        # 1b5f8 <_ZTVN3MPI4InfoE+0x10>
    be2b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    be2f:	48 89 10             	mov    %rdx,(%rax)
    be32:	90                   	nop
    be33:	5d                   	pop    %rbp
    be34:	c3                   	ret    
    be35:	90                   	nop

000000000000be36 <_ZN3MPI4InfoD0Ev>:
    be36:	f3 0f 1e fa          	endbr64 
    be3a:	55                   	push   %rbp
    be3b:	48 89 e5             	mov    %rsp,%rbp
    be3e:	48 83 ec 10          	sub    $0x10,%rsp
    be42:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    be46:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    be4a:	48 89 c7             	mov    %rax,%rdi
    be4d:	e8 c6 ff ff ff       	call   be18 <_ZN3MPI4InfoD1Ev>
    be52:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    be56:	be 10 00 00 00       	mov    $0x10,%esi
    be5b:	48 89 c7             	mov    %rax,%rdi
    be5e:	e8 9d e5 ff ff       	call   a400 <_ZdlPvm@plt>
    be63:	c9                   	leave  
    be64:	c3                   	ret    
    be65:	90                   	nop

000000000000be66 <_ZNK3MPI4InfocvP11ompi_info_tEv>:
  // comparison, don't need for info

  // inter-language operability
  Info& operator= (const MPI_Info &i) {
    mpi_info = i; return *this; }
  operator MPI_Info () const { return mpi_info; }
    be66:	f3 0f 1e fa          	endbr64 
    be6a:	55                   	push   %rbp
    be6b:	48 89 e5             	mov    %rsp,%rbp
    be6e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    be72:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    be76:	48 8b 40 08          	mov    0x8(%rax),%rax
    be7a:	5d                   	pop    %rbp
    be7b:	c3                   	ret    

000000000000be7c <_ZNK3MPI8Datatype17Create_contiguousEi>:
// Point-to-Point Communication
//

inline MPI::Datatype
MPI::Datatype::Create_contiguous(int count) const
{
    be7c:	f3 0f 1e fa          	endbr64 
    be80:	55                   	push   %rbp
    be81:	48 89 e5             	mov    %rsp,%rbp
    be84:	48 83 ec 30          	sub    $0x30,%rsp
    be88:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    be8c:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    be90:	89 55 dc             	mov    %edx,-0x24(%rbp)
    be93:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    be9a:	00 00 
    be9c:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    bea0:	31 c0                	xor    %eax,%eax
  MPI_Datatype newtype;
  (void)MPI_Type_contiguous(count, mpi_datatype, &newtype);
    bea2:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    bea6:	48 8b 48 08          	mov    0x8(%rax),%rcx
    beaa:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    beae:	8b 45 dc             	mov    -0x24(%rbp),%eax
    beb1:	48 89 ce             	mov    %rcx,%rsi
    beb4:	89 c7                	mov    %eax,%edi
    beb6:	e8 75 e2 ff ff       	call   a130 <MPI_Type_contiguous@plt>
  return newtype;
    bebb:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bebf:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    bec3:	48 89 d6             	mov    %rdx,%rsi
    bec6:	48 89 c7             	mov    %rax,%rdi
    bec9:	e8 5c f9 ff ff       	call   b82a <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
}
    bece:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bed2:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    bed9:	00 00 
    bedb:	74 05                	je     bee2 <_ZNK3MPI8Datatype17Create_contiguousEi+0x66>
    bedd:	e8 fe e5 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    bee2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    bee6:	c9                   	leave  
    bee7:	c3                   	ret    

000000000000bee8 <_ZNK3MPI8Datatype13Create_vectorEiii>:

inline MPI::Datatype
MPI::Datatype::Create_vector(int count, int blocklength,
			     int stride) const
{
    bee8:	f3 0f 1e fa          	endbr64 
    beec:	55                   	push   %rbp
    beed:	48 89 e5             	mov    %rsp,%rbp
    bef0:	48 83 ec 30          	sub    $0x30,%rsp
    bef4:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    bef8:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    befc:	89 55 dc             	mov    %edx,-0x24(%rbp)
    beff:	89 4d d8             	mov    %ecx,-0x28(%rbp)
    bf02:	44 89 45 d4          	mov    %r8d,-0x2c(%rbp)
    bf06:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    bf0d:	00 00 
    bf0f:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    bf13:	31 c0                	xor    %eax,%eax
  MPI_Datatype newtype;
  (void)MPI_Type_vector(count, blocklength, stride, mpi_datatype, &newtype);
    bf15:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    bf19:	48 8b 48 08          	mov    0x8(%rax),%rcx
    bf1d:	48 8d 7d f0          	lea    -0x10(%rbp),%rdi
    bf21:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    bf24:	8b 75 d8             	mov    -0x28(%rbp),%esi
    bf27:	8b 45 dc             	mov    -0x24(%rbp),%eax
    bf2a:	49 89 f8             	mov    %rdi,%r8
    bf2d:	89 c7                	mov    %eax,%edi
    bf2f:	e8 dc e9 ff ff       	call   a910 <MPI_Type_vector@plt>
  return newtype;
    bf34:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bf38:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    bf3c:	48 89 d6             	mov    %rdx,%rsi
    bf3f:	48 89 c7             	mov    %rax,%rdi
    bf42:	e8 e3 f8 ff ff       	call   b82a <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
}
    bf47:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bf4b:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    bf52:	00 00 
    bf54:	74 05                	je     bf5b <_ZNK3MPI8Datatype13Create_vectorEiii+0x73>
    bf56:	e8 85 e5 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    bf5b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    bf5f:	c9                   	leave  
    bf60:	c3                   	ret    
    bf61:	90                   	nop

000000000000bf62 <_ZNK3MPI8Datatype14Create_indexedEiPKiS2_>:

inline MPI::Datatype
MPI::Datatype::Create_indexed(int count,
				     const int array_of_blocklengths[],
				     const int array_of_displacements[]) const
{
    bf62:	f3 0f 1e fa          	endbr64 
    bf66:	55                   	push   %rbp
    bf67:	48 89 e5             	mov    %rsp,%rbp
    bf6a:	48 83 ec 40          	sub    $0x40,%rsp
    bf6e:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    bf72:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    bf76:	89 55 dc             	mov    %edx,-0x24(%rbp)
    bf79:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    bf7d:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    bf81:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    bf88:	00 00 
    bf8a:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    bf8e:	31 c0                	xor    %eax,%eax
  MPI_Datatype newtype;
  (void)MPI_Type_indexed(count, const_cast<int *>(array_of_blocklengths),
			 const_cast<int *>(array_of_displacements), mpi_datatype, &newtype);
    bf90:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    bf94:	48 8b 48 08          	mov    0x8(%rax),%rcx
  (void)MPI_Type_indexed(count, const_cast<int *>(array_of_blocklengths),
    bf98:	48 8d 7d f0          	lea    -0x10(%rbp),%rdi
    bf9c:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    bfa0:	48 8b 75 d0          	mov    -0x30(%rbp),%rsi
    bfa4:	8b 45 dc             	mov    -0x24(%rbp),%eax
    bfa7:	49 89 f8             	mov    %rdi,%r8
    bfaa:	89 c7                	mov    %eax,%edi
    bfac:	e8 8f e9 ff ff       	call   a940 <MPI_Type_indexed@plt>
  return newtype;
    bfb1:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bfb5:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    bfb9:	48 89 d6             	mov    %rdx,%rsi
    bfbc:	48 89 c7             	mov    %rax,%rdi
    bfbf:	e8 66 f8 ff ff       	call   b82a <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
}
    bfc4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bfc8:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    bfcf:	00 00 
    bfd1:	74 05                	je     bfd8 <_ZNK3MPI8Datatype14Create_indexedEiPKiS2_+0x76>
    bfd3:	e8 08 e5 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    bfd8:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    bfdc:	c9                   	leave  
    bfdd:	c3                   	ret    

000000000000bfde <_ZNK3MPI8Datatype15Create_hindexedEiPKiPKl>:
}

inline MPI::Datatype
MPI::Datatype::Create_hindexed(int count, const int array_of_blocklengths[],
				      const MPI::Aint array_of_displacements[]) const
{
    bfde:	f3 0f 1e fa          	endbr64 
    bfe2:	55                   	push   %rbp
    bfe3:	48 89 e5             	mov    %rsp,%rbp
    bfe6:	48 83 ec 40          	sub    $0x40,%rsp
    bfea:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    bfee:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    bff2:	89 55 dc             	mov    %edx,-0x24(%rbp)
    bff5:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    bff9:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    bffd:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c004:	00 00 
    c006:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c00a:	31 c0                	xor    %eax,%eax
  MPI_Datatype newtype;
  (void)MPI_Type_create_hindexed(count, const_cast<int *>(array_of_blocklengths),
                                 const_cast<MPI_Aint*>(array_of_displacements),
                                 mpi_datatype, &newtype) ;
    c00c:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c010:	48 8b 48 08          	mov    0x8(%rax),%rcx
  (void)MPI_Type_create_hindexed(count, const_cast<int *>(array_of_blocklengths),
    c014:	48 8d 7d f0          	lea    -0x10(%rbp),%rdi
    c018:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    c01c:	48 8b 75 d0          	mov    -0x30(%rbp),%rsi
    c020:	8b 45 dc             	mov    -0x24(%rbp),%eax
    c023:	49 89 f8             	mov    %rdi,%r8
    c026:	89 c7                	mov    %eax,%edi
    c028:	e8 e3 e4 ff ff       	call   a510 <MPI_Type_create_hindexed@plt>
  return newtype;
    c02d:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c031:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c035:	48 89 d6             	mov    %rdx,%rsi
    c038:	48 89 c7             	mov    %rax,%rdi
    c03b:	e8 ea f7 ff ff       	call   b82a <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
}
    c040:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c044:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c04b:	00 00 
    c04d:	74 05                	je     c054 <_ZNK3MPI8Datatype15Create_hindexedEiPKiPKl+0x76>
    c04f:	e8 8c e4 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    c054:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c058:	c9                   	leave  
    c059:	c3                   	ret    

000000000000c05a <_ZNK3MPI8Datatype14Create_hvectorEiil>:

inline MPI::Datatype
MPI::Datatype::Create_hvector(int count, int blocklength,
				     MPI::Aint stride) const
{
    c05a:	f3 0f 1e fa          	endbr64 
    c05e:	55                   	push   %rbp
    c05f:	48 89 e5             	mov    %rsp,%rbp
    c062:	48 83 ec 30          	sub    $0x30,%rsp
    c066:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c06a:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c06e:	89 55 dc             	mov    %edx,-0x24(%rbp)
    c071:	89 4d d8             	mov    %ecx,-0x28(%rbp)
    c074:	4c 89 45 d0          	mov    %r8,-0x30(%rbp)
    c078:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c07f:	00 00 
    c081:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c085:	31 c0                	xor    %eax,%eax
  MPI_Datatype newtype;
  (void)MPI_Type_create_hvector(count, blocklength, (MPI_Aint)stride,
                                mpi_datatype, &newtype);
    c087:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c08b:	48 8b 48 08          	mov    0x8(%rax),%rcx
  (void)MPI_Type_create_hvector(count, blocklength, (MPI_Aint)stride,
    c08f:	48 8d 7d f0          	lea    -0x10(%rbp),%rdi
    c093:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    c097:	8b 75 d8             	mov    -0x28(%rbp),%esi
    c09a:	8b 45 dc             	mov    -0x24(%rbp),%eax
    c09d:	49 89 f8             	mov    %rdi,%r8
    c0a0:	89 c7                	mov    %eax,%edi
    c0a2:	e8 b9 de ff ff       	call   9f60 <MPI_Type_create_hvector@plt>

  return newtype;
    c0a7:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c0ab:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c0af:	48 89 d6             	mov    %rdx,%rsi
    c0b2:	48 89 c7             	mov    %rax,%rdi
    c0b5:	e8 70 f7 ff ff       	call   b82a <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
}
    c0ba:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c0be:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c0c5:	00 00 
    c0c7:	74 05                	je     c0ce <_ZNK3MPI8Datatype14Create_hvectorEiil+0x74>
    c0c9:	e8 12 e4 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    c0ce:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c0d2:	c9                   	leave  
    c0d3:	c3                   	ret    

000000000000c0d4 <_ZNK3MPI8Datatype20Create_indexed_blockEiiPKi>:

inline MPI::Datatype
MPI::Datatype::Create_indexed_block(int count, int blocklength,
				    const int array_of_displacements[]) const
{
    c0d4:	f3 0f 1e fa          	endbr64 
    c0d8:	55                   	push   %rbp
    c0d9:	48 89 e5             	mov    %rsp,%rbp
    c0dc:	48 83 ec 30          	sub    $0x30,%rsp
    c0e0:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c0e4:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c0e8:	89 55 dc             	mov    %edx,-0x24(%rbp)
    c0eb:	89 4d d8             	mov    %ecx,-0x28(%rbp)
    c0ee:	4c 89 45 d0          	mov    %r8,-0x30(%rbp)
    c0f2:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c0f9:	00 00 
    c0fb:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c0ff:	31 c0                	xor    %eax,%eax
  MPI_Datatype newtype;
  (void)MPI_Type_create_indexed_block(count, blocklength, const_cast<int *>(array_of_displacements),
                                      mpi_datatype, &newtype);
    c101:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c105:	48 8b 48 08          	mov    0x8(%rax),%rcx
  (void)MPI_Type_create_indexed_block(count, blocklength, const_cast<int *>(array_of_displacements),
    c109:	48 8d 7d f0          	lea    -0x10(%rbp),%rdi
    c10d:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    c111:	8b 75 d8             	mov    -0x28(%rbp),%esi
    c114:	8b 45 dc             	mov    -0x24(%rbp),%eax
    c117:	49 89 f8             	mov    %rdi,%r8
    c11a:	89 c7                	mov    %eax,%edi
    c11c:	e8 9f e5 ff ff       	call   a6c0 <MPI_Type_create_indexed_block@plt>
  return newtype;
    c121:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c125:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c129:	48 89 d6             	mov    %rdx,%rsi
    c12c:	48 89 c7             	mov    %rax,%rdi
    c12f:	e8 f6 f6 ff ff       	call   b82a <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
}
    c134:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c138:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c13f:	00 00 
    c141:	74 05                	je     c148 <_ZNK3MPI8Datatype20Create_indexed_blockEiiPKi+0x74>
    c143:	e8 98 e3 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    c148:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c14c:	c9                   	leave  
    c14d:	c3                   	ret    

000000000000c14e <_ZNK3MPI8Datatype14Create_resizedEll>:

inline MPI::Datatype
MPI::Datatype::Create_resized(const MPI::Aint lb, const MPI::Aint extent) const
{
    c14e:	f3 0f 1e fa          	endbr64 
    c152:	55                   	push   %rbp
    c153:	48 89 e5             	mov    %rsp,%rbp
    c156:	48 83 ec 30          	sub    $0x30,%rsp
    c15a:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c15e:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c162:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    c166:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    c16a:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c171:	00 00 
    c173:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c177:	31 c0                	xor    %eax,%eax
    MPI_Datatype newtype;

    (void) MPI_Type_create_resized(mpi_datatype, lb, extent, &newtype);
    c179:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c17d:	48 8b 40 08          	mov    0x8(%rax),%rax
    c181:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    c185:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    c189:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    c18d:	48 89 c7             	mov    %rax,%rdi
    c190:	e8 0b e8 ff ff       	call   a9a0 <MPI_Type_create_resized@plt>
    return(newtype);
    c195:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c199:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c19d:	48 89 d6             	mov    %rdx,%rsi
    c1a0:	48 89 c7             	mov    %rax,%rdi
    c1a3:	e8 82 f6 ff ff       	call   b82a <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
}
    c1a8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c1ac:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c1b3:	00 00 
    c1b5:	74 05                	je     c1bc <_ZNK3MPI8Datatype14Create_resizedEll+0x6e>
    c1b7:	e8 24 e3 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    c1bc:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c1c0:	c9                   	leave  
    c1c1:	c3                   	ret    

000000000000c1c2 <_ZNK3MPI8Datatype8Get_sizeEv>:

inline int
MPI::Datatype::Get_size() const
{
    c1c2:	f3 0f 1e fa          	endbr64 
    c1c6:	55                   	push   %rbp
    c1c7:	48 89 e5             	mov    %rsp,%rbp
    c1ca:	48 83 ec 20          	sub    $0x20,%rsp
    c1ce:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c1d2:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c1d9:	00 00 
    c1db:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c1df:	31 c0                	xor    %eax,%eax
  int size;
  (void)MPI_Type_size(mpi_datatype, &size);
    c1e1:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c1e5:	48 8b 40 08          	mov    0x8(%rax),%rax
    c1e9:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    c1ed:	48 89 d6             	mov    %rdx,%rsi
    c1f0:	48 89 c7             	mov    %rax,%rdi
    c1f3:	e8 e8 e3 ff ff       	call   a5e0 <MPI_Type_size@plt>
  return size;
    c1f8:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    c1fb:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    c1ff:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    c206:	00 00 
    c208:	74 05                	je     c20f <_ZNK3MPI8Datatype8Get_sizeEv+0x4d>
    c20a:	e8 d1 e2 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    c20f:	c9                   	leave  
    c210:	c3                   	ret    
    c211:	90                   	nop

000000000000c212 <_ZNK3MPI8Datatype10Get_extentERlS1_>:

inline void
MPI::Datatype::Get_extent(MPI::Aint& lb, MPI::Aint& extent) const
{
    c212:	f3 0f 1e fa          	endbr64 
    c216:	55                   	push   %rbp
    c217:	48 89 e5             	mov    %rsp,%rbp
    c21a:	48 83 ec 20          	sub    $0x20,%rsp
    c21e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c222:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c226:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
  (void)MPI_Type_get_extent(mpi_datatype, &lb, &extent);
    c22a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c22e:	48 8b 40 08          	mov    0x8(%rax),%rax
    c232:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    c236:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    c23a:	48 89 ce             	mov    %rcx,%rsi
    c23d:	48 89 c7             	mov    %rax,%rdi
    c240:	e8 ab dd ff ff       	call   9ff0 <MPI_Type_get_extent@plt>
}
    c245:	90                   	nop
    c246:	c9                   	leave  
    c247:	c3                   	ret    

000000000000c248 <_ZNK3MPI8Datatype15Get_true_extentERlS1_>:

inline void
MPI::Datatype::Get_true_extent(MPI::Aint& lb, MPI::Aint& extent) const
{
    c248:	f3 0f 1e fa          	endbr64 
    c24c:	55                   	push   %rbp
    c24d:	48 89 e5             	mov    %rsp,%rbp
    c250:	48 83 ec 20          	sub    $0x20,%rsp
    c254:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c258:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c25c:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    (void) MPI_Type_get_true_extent(mpi_datatype, &lb, &extent);
    c260:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c264:	48 8b 40 08          	mov    0x8(%rax),%rax
    c268:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    c26c:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    c270:	48 89 ce             	mov    %rcx,%rsi
    c273:	48 89 c7             	mov    %rax,%rdi
    c276:	e8 85 de ff ff       	call   a100 <MPI_Type_get_true_extent@plt>
}
    c27b:	90                   	nop
    c27c:	c9                   	leave  
    c27d:	c3                   	ret    

000000000000c27e <_ZN3MPI8Datatype6CommitEv>:

inline void
MPI::Datatype::Commit()
{
    c27e:	f3 0f 1e fa          	endbr64 
    c282:	55                   	push   %rbp
    c283:	48 89 e5             	mov    %rsp,%rbp
    c286:	48 83 ec 10          	sub    $0x10,%rsp
    c28a:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
  (void)MPI_Type_commit(&mpi_datatype);
    c28e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c292:	48 83 c0 08          	add    $0x8,%rax
    c296:	48 89 c7             	mov    %rax,%rdi
    c299:	e8 12 e3 ff ff       	call   a5b0 <MPI_Type_commit@plt>
}
    c29e:	90                   	nop
    c29f:	c9                   	leave  
    c2a0:	c3                   	ret    
    c2a1:	90                   	nop

000000000000c2a2 <_ZNK3MPI8Datatype4PackEPKviPviRiRKNS_4CommE>:

inline void
MPI::Datatype::Pack(const void* inbuf, int incount,
			   void *outbuf, int outsize,
			   int& position, const MPI::Comm &comm) const
{
    c2a2:	f3 0f 1e fa          	endbr64 
    c2a6:	55                   	push   %rbp
    c2a7:	48 89 e5             	mov    %rsp,%rbp
    c2aa:	48 83 ec 30          	sub    $0x30,%rsp
    c2ae:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c2b2:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c2b6:	89 55 ec             	mov    %edx,-0x14(%rbp)
    c2b9:	48 89 4d e0          	mov    %rcx,-0x20(%rbp)
    c2bd:	44 89 45 e8          	mov    %r8d,-0x18(%rbp)
    c2c1:	4c 89 4d d8          	mov    %r9,-0x28(%rbp)
  (void)MPI_Pack(const_cast<void *>(inbuf), incount,  mpi_datatype, outbuf,
    c2c5:	48 8b 45 10          	mov    0x10(%rbp),%rax
    c2c9:	48 89 c7             	mov    %rax,%rdi
    c2cc:	e8 c9 f8 ff ff       	call   bb9a <_ZNK3MPI9Comm_NullcvP19ompi_communicator_tEv>
    c2d1:	48 89 c7             	mov    %rax,%rdi
    c2d4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c2d8:	48 8b 50 08          	mov    0x8(%rax),%rdx
    c2dc:	4c 8b 4d d8          	mov    -0x28(%rbp),%r9
    c2e0:	44 8b 45 e8          	mov    -0x18(%rbp),%r8d
    c2e4:	48 8b 4d e0          	mov    -0x20(%rbp),%rcx
    c2e8:	8b 75 ec             	mov    -0x14(%rbp),%esi
    c2eb:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    c2ef:	48 83 ec 08          	sub    $0x8,%rsp
    c2f3:	57                   	push   %rdi
    c2f4:	48 89 c7             	mov    %rax,%rdi
    c2f7:	e8 b4 e0 ff ff       	call   a3b0 <MPI_Pack@plt>
    c2fc:	48 83 c4 10          	add    $0x10,%rsp
		 outsize, &position, comm);
}
    c300:	90                   	nop
    c301:	c9                   	leave  
    c302:	c3                   	ret    
    c303:	90                   	nop

000000000000c304 <_ZNK3MPI8Datatype6UnpackEPKviPviRiRKNS_4CommE>:

inline void
MPI::Datatype::Unpack(const void* inbuf, int insize,
			     void *outbuf, int outcount, int& position,
			     const MPI::Comm& comm) const
{
    c304:	f3 0f 1e fa          	endbr64 
    c308:	55                   	push   %rbp
    c309:	48 89 e5             	mov    %rsp,%rbp
    c30c:	48 83 ec 30          	sub    $0x30,%rsp
    c310:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c314:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c318:	89 55 ec             	mov    %edx,-0x14(%rbp)
    c31b:	48 89 4d e0          	mov    %rcx,-0x20(%rbp)
    c31f:	44 89 45 e8          	mov    %r8d,-0x18(%rbp)
    c323:	4c 89 4d d8          	mov    %r9,-0x28(%rbp)
  (void)MPI_Unpack(const_cast<void *>(inbuf), insize, &position,
    c327:	48 8b 45 10          	mov    0x10(%rbp),%rax
    c32b:	48 89 c7             	mov    %rax,%rdi
    c32e:	e8 67 f8 ff ff       	call   bb9a <_ZNK3MPI9Comm_NullcvP19ompi_communicator_tEv>
    c333:	48 89 c7             	mov    %rax,%rdi
		   outbuf, outcount, mpi_datatype, comm);
    c336:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c33a:	4c 8b 48 08          	mov    0x8(%rax),%r9
  (void)MPI_Unpack(const_cast<void *>(inbuf), insize, &position,
    c33e:	44 8b 45 e8          	mov    -0x18(%rbp),%r8d
    c342:	48 8b 4d e0          	mov    -0x20(%rbp),%rcx
    c346:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
    c34a:	8b 75 ec             	mov    -0x14(%rbp),%esi
    c34d:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    c351:	48 83 ec 08          	sub    $0x8,%rsp
    c355:	57                   	push   %rdi
    c356:	48 89 c7             	mov    %rax,%rdi
    c359:	e8 c2 dc ff ff       	call   a020 <MPI_Unpack@plt>
    c35e:	48 83 c4 10          	add    $0x10,%rsp
}
    c362:	90                   	nop
    c363:	c9                   	leave  
    c364:	c3                   	ret    
    c365:	90                   	nop

000000000000c366 <_ZNK3MPI8Datatype9Pack_sizeEiRKNS_4CommE>:

inline int
MPI::Datatype::Pack_size(int incount, const MPI::Comm& comm) const
{
    c366:	f3 0f 1e fa          	endbr64 
    c36a:	55                   	push   %rbp
    c36b:	48 89 e5             	mov    %rsp,%rbp
    c36e:	48 83 ec 30          	sub    $0x30,%rsp
    c372:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c376:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    c379:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    c37d:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c384:	00 00 
    c386:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c38a:	31 c0                	xor    %eax,%eax
  int size;
  (void)MPI_Pack_size(incount, mpi_datatype, comm, &size);
    c38c:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    c390:	48 89 c7             	mov    %rax,%rdi
    c393:	e8 02 f8 ff ff       	call   bb9a <_ZNK3MPI9Comm_NullcvP19ompi_communicator_tEv>
    c398:	48 89 c7             	mov    %rax,%rdi
    c39b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c39f:	48 8b 70 08          	mov    0x8(%rax),%rsi
    c3a3:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    c3a7:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    c3aa:	48 89 d1             	mov    %rdx,%rcx
    c3ad:	48 89 fa             	mov    %rdi,%rdx
    c3b0:	89 c7                	mov    %eax,%edi
    c3b2:	e8 79 db ff ff       	call   9f30 <MPI_Pack_size@plt>
  return size;
    c3b7:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    c3ba:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    c3be:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    c3c5:	00 00 
    c3c7:	74 05                	je     c3ce <_ZNK3MPI8Datatype9Pack_sizeEiRKNS_4CommE+0x68>
    c3c9:	e8 12 e1 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    c3ce:	c9                   	leave  
    c3cf:	c3                   	ret    

000000000000c3d0 <_ZNK3MPI8Datatype13Pack_externalEPKcPKviPvlRl>:

inline void
MPI::Datatype::Pack_external(const char* datarep, const void* inbuf, int incount,
            void* outbuf, MPI::Aint outsize, MPI::Aint& position) const
{
    c3d0:	f3 0f 1e fa          	endbr64 
    c3d4:	55                   	push   %rbp
    c3d5:	48 89 e5             	mov    %rsp,%rbp
    c3d8:	48 83 ec 30          	sub    $0x30,%rsp
    c3dc:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c3e0:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c3e4:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    c3e8:	89 4d e4             	mov    %ecx,-0x1c(%rbp)
    c3eb:	4c 89 45 d8          	mov    %r8,-0x28(%rbp)
    c3ef:	4c 89 4d d0          	mov    %r9,-0x30(%rbp)
    (void)MPI_Pack_external(const_cast<char *>(datarep), const_cast<void *>(inbuf),
                             incount, mpi_datatype, outbuf, outsize, &position);
    c3f3:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c3f7:	48 8b 48 08          	mov    0x8(%rax),%rcx
    (void)MPI_Pack_external(const_cast<char *>(datarep), const_cast<void *>(inbuf),
    c3fb:	4c 8b 45 d0          	mov    -0x30(%rbp),%r8
    c3ff:	48 8b 7d d8          	mov    -0x28(%rbp),%rdi
    c403:	8b 55 e4             	mov    -0x1c(%rbp),%edx
    c406:	48 8b 75 e8          	mov    -0x18(%rbp),%rsi
    c40a:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    c40e:	48 83 ec 08          	sub    $0x8,%rsp
    c412:	ff 75 10             	push   0x10(%rbp)
    c415:	4d 89 c1             	mov    %r8,%r9
    c418:	49 89 f8             	mov    %rdi,%r8
    c41b:	48 89 c7             	mov    %rax,%rdi
    c41e:	e8 dd dd ff ff       	call   a200 <MPI_Pack_external@plt>
    c423:	48 83 c4 10          	add    $0x10,%rsp
}
    c427:	90                   	nop
    c428:	c9                   	leave  
    c429:	c3                   	ret    

000000000000c42a <_ZNK3MPI8Datatype18Pack_external_sizeEPKci>:

inline MPI::Aint
MPI::Datatype::Pack_external_size(const char* datarep, int incount) const
{
    c42a:	f3 0f 1e fa          	endbr64 
    c42e:	55                   	push   %rbp
    c42f:	48 89 e5             	mov    %rsp,%rbp
    c432:	48 83 ec 30          	sub    $0x30,%rsp
    c436:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c43a:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c43e:	89 55 dc             	mov    %edx,-0x24(%rbp)
    c441:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c448:	00 00 
    c44a:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c44e:	31 c0                	xor    %eax,%eax
    MPI_Aint addr;
    (void)MPI_Pack_external_size(const_cast<char *>(datarep), incount, mpi_datatype, &addr);
    c450:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c454:	48 8b 50 08          	mov    0x8(%rax),%rdx
    c458:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    c45c:	8b 75 dc             	mov    -0x24(%rbp),%esi
    c45f:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c463:	48 89 c7             	mov    %rax,%rdi
    c466:	e8 65 e6 ff ff       	call   aad0 <MPI_Pack_external_size@plt>
    return addr;
    c46b:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
}
    c46f:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    c473:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    c47a:	00 00 
    c47c:	74 05                	je     c483 <_ZNK3MPI8Datatype18Pack_external_sizeEPKci+0x59>
    c47e:	e8 5d e0 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    c483:	c9                   	leave  
    c484:	c3                   	ret    
    c485:	90                   	nop

000000000000c486 <_ZNK3MPI8Datatype15Unpack_externalEPKcPKvlRlPvi>:

inline void
MPI::Datatype::Unpack_external(const char* datarep, const void* inbuf,
            MPI::Aint insize, MPI::Aint& position, void* outbuf, int outcount) const
{
    c486:	f3 0f 1e fa          	endbr64 
    c48a:	55                   	push   %rbp
    c48b:	48 89 e5             	mov    %rsp,%rbp
    c48e:	48 83 ec 30          	sub    $0x30,%rsp
    c492:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c496:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c49a:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    c49e:	48 89 4d e0          	mov    %rcx,-0x20(%rbp)
    c4a2:	4c 89 45 d8          	mov    %r8,-0x28(%rbp)
    c4a6:	4c 89 4d d0          	mov    %r9,-0x30(%rbp)
    (void)MPI_Unpack_external(const_cast<char *>(datarep), const_cast<void *>(inbuf),
                               insize, &position, outbuf, outcount, mpi_datatype);
    c4aa:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c4ae:	48 8b 78 08          	mov    0x8(%rax),%rdi
    (void)MPI_Unpack_external(const_cast<char *>(datarep), const_cast<void *>(inbuf),
    c4b2:	4c 8b 45 d0          	mov    -0x30(%rbp),%r8
    c4b6:	48 8b 4d d8          	mov    -0x28(%rbp),%rcx
    c4ba:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    c4be:	48 8b 75 e8          	mov    -0x18(%rbp),%rsi
    c4c2:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    c4c6:	48 83 ec 08          	sub    $0x8,%rsp
    c4ca:	57                   	push   %rdi
    c4cb:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
    c4cf:	48 89 c7             	mov    %rax,%rdi
    c4d2:	e8 99 da ff ff       	call   9f70 <MPI_Unpack_external@plt>
    c4d7:	48 83 c4 10          	add    $0x10,%rsp
}
    c4db:	90                   	nop
    c4dc:	c9                   	leave  
    c4dd:	c3                   	ret    

000000000000c4de <_ZNK3MPI8Datatype15Create_subarrayEiPKiS2_S2_i>:
inline MPI::Datatype
MPI::Datatype::Create_subarray(int ndims, const int array_of_sizes[],
				      const int array_of_subsizes[],
				      const int array_of_starts[], int order)
  const
{
    c4de:	f3 0f 1e fa          	endbr64 
    c4e2:	55                   	push   %rbp
    c4e3:	48 89 e5             	mov    %rsp,%rbp
    c4e6:	48 83 ec 40          	sub    $0x40,%rsp
    c4ea:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c4ee:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c4f2:	89 55 dc             	mov    %edx,-0x24(%rbp)
    c4f5:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    c4f9:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    c4fd:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    c501:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c508:	00 00 
    c50a:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c50e:	31 c0                	xor    %eax,%eax
  MPI_Datatype type;
  (void) MPI_Type_create_subarray(ndims, const_cast<int *>(array_of_sizes),
                                  const_cast<int *>(array_of_subsizes),
                                  const_cast<int *>(array_of_starts),
                                  order, mpi_datatype, &type);
    c510:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c514:	4c 8b 40 08          	mov    0x8(%rax),%r8
  (void) MPI_Type_create_subarray(ndims, const_cast<int *>(array_of_sizes),
    c518:	48 8b 4d c0          	mov    -0x40(%rbp),%rcx
    c51c:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    c520:	48 8b 75 d0          	mov    -0x30(%rbp),%rsi
    c524:	8b 45 dc             	mov    -0x24(%rbp),%eax
    c527:	48 83 ec 08          	sub    $0x8,%rsp
    c52b:	48 8d 7d f0          	lea    -0x10(%rbp),%rdi
    c52f:	57                   	push   %rdi
    c530:	4d 89 c1             	mov    %r8,%r9
    c533:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    c537:	89 c7                	mov    %eax,%edi
    c539:	e8 72 e3 ff ff       	call   a8b0 <MPI_Type_create_subarray@plt>
    c53e:	48 83 c4 10          	add    $0x10,%rsp
  return type;
    c542:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c546:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c54a:	48 89 d6             	mov    %rdx,%rsi
    c54d:	48 89 c7             	mov    %rax,%rdi
    c550:	e8 d5 f2 ff ff       	call   b82a <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
}
    c555:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c559:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c560:	00 00 
    c562:	74 05                	je     c569 <_ZNK3MPI8Datatype15Create_subarrayEiPKiS2_S2_i+0x8b>
    c564:	e8 77 df ff ff       	call   a4e0 <__stack_chk_fail@plt>
    c569:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c56d:	c9                   	leave  
    c56e:	c3                   	ret    
    c56f:	90                   	nop

000000000000c570 <_ZNK3MPI8Datatype13Create_darrayEiiiPKiS2_S2_S2_i>:
inline MPI::Datatype
MPI::Datatype::Create_darray(int size, int rank, int ndims,
                   const int array_of_gsizes[], const int array_of_distribs[],
                   const int array_of_dargs[],  const int array_of_psizes[],
                   int order) const
{
    c570:	f3 0f 1e fa          	endbr64 
    c574:	55                   	push   %rbp
    c575:	48 89 e5             	mov    %rsp,%rbp
    c578:	48 83 ec 50          	sub    $0x50,%rsp
    c57c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c580:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c584:	89 55 dc             	mov    %edx,-0x24(%rbp)
    c587:	89 4d d8             	mov    %ecx,-0x28(%rbp)
    c58a:	44 89 45 d4          	mov    %r8d,-0x2c(%rbp)
    c58e:	4c 89 4d c8          	mov    %r9,-0x38(%rbp)
    c592:	48 8b 45 10          	mov    0x10(%rbp),%rax
    c596:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
    c59a:	48 8b 45 18          	mov    0x18(%rbp),%rax
    c59e:	48 89 45 b8          	mov    %rax,-0x48(%rbp)
    c5a2:	48 8b 45 20          	mov    0x20(%rbp),%rax
    c5a6:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
    c5aa:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c5b1:	00 00 
    c5b3:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c5b7:	31 c0                	xor    %eax,%eax
    (void) MPI_Type_create_darray(size, rank, ndims,
                   const_cast<int *>(array_of_gsizes),
                   const_cast<int *>(array_of_distribs),
                   const_cast<int *>(array_of_dargs),
                   const_cast<int *>(array_of_psizes),
                   order, mpi_datatype, &type);
    c5b9:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c5bd:	48 8b 78 08          	mov    0x8(%rax),%rdi
    (void) MPI_Type_create_darray(size, rank, ndims,
    c5c1:	4c 8b 4d b8          	mov    -0x48(%rbp),%r9
    c5c5:	4c 8b 55 c0          	mov    -0x40(%rbp),%r10
    c5c9:	48 8b 4d c8          	mov    -0x38(%rbp),%rcx
    c5cd:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    c5d0:	8b 75 d8             	mov    -0x28(%rbp),%esi
    c5d3:	8b 45 dc             	mov    -0x24(%rbp),%eax
    c5d6:	4c 8d 45 f0          	lea    -0x10(%rbp),%r8
    c5da:	41 50                	push   %r8
    c5dc:	57                   	push   %rdi
    c5dd:	8b 7d 28             	mov    0x28(%rbp),%edi
    c5e0:	57                   	push   %rdi
    c5e1:	ff 75 b0             	push   -0x50(%rbp)
    c5e4:	4d 89 d0             	mov    %r10,%r8
    c5e7:	89 c7                	mov    %eax,%edi
    c5e9:	e8 c2 dc ff ff       	call   a2b0 <MPI_Type_create_darray@plt>
    c5ee:	48 83 c4 20          	add    $0x20,%rsp
    return type;
    c5f2:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c5f6:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c5fa:	48 89 d6             	mov    %rdx,%rsi
    c5fd:	48 89 c7             	mov    %rax,%rdi
    c600:	e8 25 f2 ff ff       	call   b82a <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
}
    c605:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c609:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c610:	00 00 
    c612:	74 05                	je     c619 <_ZNK3MPI8Datatype13Create_darrayEiiiPKiS2_S2_S2_i+0xa9>
    c614:	e8 c7 de ff ff       	call   a4e0 <__stack_chk_fail@plt>
    c619:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c61d:	c9                   	leave  
    c61e:	c3                   	ret    
    c61f:	90                   	nop

000000000000c620 <_ZNK3MPI8Datatype3DupEv>:
//


inline MPI::Datatype
MPI::Datatype::Dup() const
{
    c620:	f3 0f 1e fa          	endbr64 
    c624:	55                   	push   %rbp
    c625:	48 89 e5             	mov    %rsp,%rbp
    c628:	48 83 ec 20          	sub    $0x20,%rsp
    c62c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c630:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c634:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c63b:	00 00 
    c63d:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c641:	31 c0                	xor    %eax,%eax
  MPI_Datatype type;
  (void) MPI_Type_dup(mpi_datatype, &type);
    c643:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c647:	48 8b 40 08          	mov    0x8(%rax),%rax
    c64b:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    c64f:	48 89 d6             	mov    %rdx,%rsi
    c652:	48 89 c7             	mov    %rax,%rdi
    c655:	e8 66 e1 ff ff       	call   a7c0 <MPI_Type_dup@plt>
  return type;
    c65a:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c65e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c662:	48 89 d6             	mov    %rdx,%rsi
    c665:	48 89 c7             	mov    %rax,%rdi
    c668:	e8 bd f1 ff ff       	call   b82a <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
}
    c66d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c671:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c678:	00 00 
    c67a:	74 05                	je     c681 <_ZNK3MPI8Datatype3DupEv+0x61>
    c67c:	e8 5f de ff ff       	call   a4e0 <__stack_chk_fail@plt>
    c681:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c685:	c9                   	leave  
    c686:	c3                   	ret    
    c687:	90                   	nop

000000000000c688 <_ZN3MPI8Datatype11Delete_attrEi>:
    return (MPI_SUCCESS == ret) ? keyval : ret;
}

inline void
MPI::Datatype::Delete_attr(int type_keyval)
{
    c688:	f3 0f 1e fa          	endbr64 
    c68c:	55                   	push   %rbp
    c68d:	48 89 e5             	mov    %rsp,%rbp
    c690:	48 83 ec 10          	sub    $0x10,%rsp
    c694:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c698:	89 75 f4             	mov    %esi,-0xc(%rbp)
  (void) MPI_Type_delete_attr(mpi_datatype, type_keyval);
    c69b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c69f:	48 8b 40 08          	mov    0x8(%rax),%rax
    c6a3:	8b 55 f4             	mov    -0xc(%rbp),%edx
    c6a6:	89 d6                	mov    %edx,%esi
    c6a8:	48 89 c7             	mov    %rax,%rdi
    c6ab:	e8 60 e3 ff ff       	call   aa10 <MPI_Type_delete_attr@plt>
}
    c6b0:	90                   	nop
    c6b1:	c9                   	leave  
    c6b2:	c3                   	ret    
    c6b3:	90                   	nop

000000000000c6b4 <_ZNK3MPI8Datatype8Get_attrEiPv>:
}

inline bool
MPI::Datatype::Get_attr(int type_keyval,
                          void* attribute_val) const
{
    c6b4:	f3 0f 1e fa          	endbr64 
    c6b8:	55                   	push   %rbp
    c6b9:	48 89 e5             	mov    %rsp,%rbp
    c6bc:	48 83 ec 30          	sub    $0x30,%rsp
    c6c0:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c6c4:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    c6c7:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    c6cb:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c6d2:	00 00 
    c6d4:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c6d8:	31 c0                	xor    %eax,%eax
  int ret;
  (void) MPI_Type_get_attr(mpi_datatype, type_keyval, attribute_val, &ret);
    c6da:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c6de:	48 8b 40 08          	mov    0x8(%rax),%rax
    c6e2:	48 8d 4d f4          	lea    -0xc(%rbp),%rcx
    c6e6:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
    c6ea:	8b 75 e4             	mov    -0x1c(%rbp),%esi
    c6ed:	48 89 c7             	mov    %rax,%rdi
    c6f0:	e8 4b dc ff ff       	call   a340 <MPI_Type_get_attr@plt>
  return OPAL_INT_TO_BOOL(ret);
    c6f5:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c6f8:	85 c0                	test   %eax,%eax
    c6fa:	0f 95 c0             	setne  %al
}
    c6fd:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    c701:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    c708:	00 00 
    c70a:	74 05                	je     c711 <_ZNK3MPI8Datatype8Get_attrEiPv+0x5d>
    c70c:	e8 cf dd ff ff       	call   a4e0 <__stack_chk_fail@plt>
    c711:	c9                   	leave  
    c712:	c3                   	ret    
    c713:	90                   	nop

000000000000c714 <_ZNK3MPI8Datatype12Get_contentsEiiiPiPlPS0_>:
inline void
MPI::Datatype::Get_contents(int max_integers, int max_addresses,
                            int max_datatypes, int array_of_integers[],
                            MPI::Aint array_of_addresses[],
                            MPI::Datatype array_of_datatypes[]) const
{
    c714:	f3 0f 1e fa          	endbr64 
    c718:	55                   	push   %rbp
    c719:	48 89 e5             	mov    %rsp,%rbp
    c71c:	48 83 ec 40          	sub    $0x40,%rsp
    c720:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c724:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    c727:	89 55 e0             	mov    %edx,-0x20(%rbp)
    c72a:	89 4d dc             	mov    %ecx,-0x24(%rbp)
    c72d:	4c 89 45 d0          	mov    %r8,-0x30(%rbp)
    c731:	4c 89 4d c8          	mov    %r9,-0x38(%rbp)
    int i;
    MPI_Datatype *c_datatypes = new MPI_Datatype[max_datatypes];
    c735:	8b 45 dc             	mov    -0x24(%rbp),%eax
    c738:	48 98                	cltq   
    c73a:	48 ba ff ff ff ff ff 	movabs $0xfffffffffffffff,%rdx
    c741:	ff ff 0f 
    c744:	48 39 d0             	cmp    %rdx,%rax
    c747:	77 4b                	ja     c794 <_ZNK3MPI8Datatype12Get_contentsEiiiPiPlPS0_+0x80>
    c749:	48 c1 e0 03          	shl    $0x3,%rax
    c74d:	48 89 c7             	mov    %rax,%rdi
    c750:	e8 6b d6 ff ff       	call   9dc0 <_Znam@plt>
    c755:	48 89 45 f8          	mov    %rax,-0x8(%rbp)

    (void) MPI_Type_get_contents(mpi_datatype, max_integers, max_addresses,
    c759:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c75d:	48 8b 40 08          	mov    0x8(%rax),%rax
    c761:	4c 8b 45 c8          	mov    -0x38(%rbp),%r8
    c765:	48 8b 7d d0          	mov    -0x30(%rbp),%rdi
    c769:	8b 4d dc             	mov    -0x24(%rbp),%ecx
    c76c:	8b 55 e0             	mov    -0x20(%rbp),%edx
    c76f:	8b 75 e4             	mov    -0x1c(%rbp),%esi
    c772:	48 83 ec 08          	sub    $0x8,%rsp
    c776:	ff 75 f8             	push   -0x8(%rbp)
    c779:	4d 89 c1             	mov    %r8,%r9
    c77c:	49 89 f8             	mov    %rdi,%r8
    c77f:	48 89 c7             	mov    %rax,%rdi
    c782:	e8 49 df ff ff       	call   a6d0 <MPI_Type_get_contents@plt>
    c787:	48 83 c4 10          	add    $0x10,%rsp
                                 const_cast<int *>(array_of_integers),
                                 const_cast<MPI_Aint*>(array_of_addresses),
                                 c_datatypes);
    // Convert the C MPI_Datatypes to the user's OUT MPI::Datatype
    // array parameter
    for (i = 0; i < max_datatypes; ++i) {
    c78b:	c7 45 f4 00 00 00 00 	movl   $0x0,-0xc(%rbp)
    c792:	eb 3b                	jmp    c7cf <_ZNK3MPI8Datatype12Get_contentsEiiiPiPlPS0_+0xbb>
    MPI_Datatype *c_datatypes = new MPI_Datatype[max_datatypes];
    c794:	e8 57 df ff ff       	call   a6f0 <__cxa_throw_bad_array_new_length@plt>
        array_of_datatypes[i] = c_datatypes[i];
    c799:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c79c:	48 98                	cltq   
    c79e:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
    c7a5:	00 
    c7a6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c7aa:	48 01 c2             	add    %rax,%rdx
    c7ad:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c7b0:	48 98                	cltq   
    c7b2:	48 c1 e0 04          	shl    $0x4,%rax
    c7b6:	48 89 c1             	mov    %rax,%rcx
    c7b9:	48 8b 45 10          	mov    0x10(%rbp),%rax
    c7bd:	48 01 c8             	add    %rcx,%rax
    c7c0:	48 89 d6             	mov    %rdx,%rsi
    c7c3:	48 89 c7             	mov    %rax,%rdi
    c7c6:	e8 8d f0 ff ff       	call   b858 <_ZN3MPI8DatatypeaSERKP15ompi_datatype_t>
    for (i = 0; i < max_datatypes; ++i) {
    c7cb:	83 45 f4 01          	addl   $0x1,-0xc(%rbp)
    c7cf:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c7d2:	3b 45 dc             	cmp    -0x24(%rbp),%eax
    c7d5:	7c c2                	jl     c799 <_ZNK3MPI8Datatype12Get_contentsEiiiPiPlPS0_+0x85>
    }
    delete[] c_datatypes;
    c7d7:	48 83 7d f8 00       	cmpq   $0x0,-0x8(%rbp)
    c7dc:	74 0c                	je     c7ea <_ZNK3MPI8Datatype12Get_contentsEiiiPiPlPS0_+0xd6>
    c7de:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c7e2:	48 89 c7             	mov    %rax,%rdi
    c7e5:	e8 e6 dd ff ff       	call   a5d0 <_ZdaPv@plt>
}
    c7ea:	90                   	nop
    c7eb:	c9                   	leave  
    c7ec:	c3                   	ret    
    c7ed:	90                   	nop

000000000000c7ee <_ZNK3MPI8Datatype12Get_envelopeERiS1_S1_S1_>:

inline void
MPI::Datatype::Get_envelope(int& num_integers, int& num_addresses,
			  int& num_datatypes, int& combiner) const
{
    c7ee:	f3 0f 1e fa          	endbr64 
    c7f2:	55                   	push   %rbp
    c7f3:	48 89 e5             	mov    %rsp,%rbp
    c7f6:	48 83 ec 30          	sub    $0x30,%rsp
    c7fa:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c7fe:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c802:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    c806:	48 89 4d e0          	mov    %rcx,-0x20(%rbp)
    c80a:	4c 89 45 d8          	mov    %r8,-0x28(%rbp)
  (void) MPI_Type_get_envelope(mpi_datatype, &num_integers, &num_addresses,
    c80e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c812:	48 8b 40 08          	mov    0x8(%rax),%rax
    c816:	48 8b 7d d8          	mov    -0x28(%rbp),%rdi
    c81a:	48 8b 4d e0          	mov    -0x20(%rbp),%rcx
    c81e:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    c822:	48 8b 75 f0          	mov    -0x10(%rbp),%rsi
    c826:	49 89 f8             	mov    %rdi,%r8
    c829:	48 89 c7             	mov    %rax,%rdi
    c82c:	e8 1f d8 ff ff       	call   a050 <MPI_Type_get_envelope@plt>
				&num_datatypes, &combiner);
}
    c831:	90                   	nop
    c832:	c9                   	leave  
    c833:	c3                   	ret    

000000000000c834 <_ZNK3MPI8Datatype8Get_nameEPcRi>:

inline void
MPI::Datatype::Get_name(char* type_name, int& resultlen) const
{
    c834:	f3 0f 1e fa          	endbr64 
    c838:	55                   	push   %rbp
    c839:	48 89 e5             	mov    %rsp,%rbp
    c83c:	48 83 ec 20          	sub    $0x20,%rsp
    c840:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c844:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c848:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
  (void) MPI_Type_get_name(mpi_datatype, type_name, &resultlen);
    c84c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c850:	48 8b 40 08          	mov    0x8(%rax),%rax
    c854:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    c858:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    c85c:	48 89 ce             	mov    %rcx,%rsi
    c85f:	48 89 c7             	mov    %rax,%rdi
    c862:	e8 49 de ff ff       	call   a6b0 <MPI_Type_get_name@plt>
}
    c867:	90                   	nop
    c868:	c9                   	leave  
    c869:	c3                   	ret    

000000000000c86a <_ZN3MPI8Datatype8Set_attrEiPKv>:

inline void
MPI::Datatype::Set_attr(int type_keyval, const void* attribute_val)
{
    c86a:	f3 0f 1e fa          	endbr64 
    c86e:	55                   	push   %rbp
    c86f:	48 89 e5             	mov    %rsp,%rbp
    c872:	48 83 ec 20          	sub    $0x20,%rsp
    c876:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c87a:	89 75 f4             	mov    %esi,-0xc(%rbp)
    c87d:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
  (void) MPI_Type_set_attr(mpi_datatype, type_keyval, const_cast<void *>(attribute_val));
    c881:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c885:	48 8b 40 08          	mov    0x8(%rax),%rax
    c889:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    c88d:	8b 4d f4             	mov    -0xc(%rbp),%ecx
    c890:	89 ce                	mov    %ecx,%esi
    c892:	48 89 c7             	mov    %rax,%rdi
    c895:	e8 f6 e1 ff ff       	call   aa90 <MPI_Type_set_attr@plt>
}
    c89a:	90                   	nop
    c89b:	c9                   	leave  
    c89c:	c3                   	ret    
    c89d:	90                   	nop

000000000000c89e <_ZN3MPI8Datatype8Set_nameEPKc>:

inline void
MPI::Datatype::Set_name(const char* type_name)
{
    c89e:	f3 0f 1e fa          	endbr64 
    c8a2:	55                   	push   %rbp
    c8a3:	48 89 e5             	mov    %rsp,%rbp
    c8a6:	48 83 ec 10          	sub    $0x10,%rsp
    c8aa:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c8ae:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
  (void) MPI_Type_set_name(mpi_datatype, const_cast<char *>(type_name));
    c8b2:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c8b6:	48 8b 40 08          	mov    0x8(%rax),%rax
    c8ba:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c8be:	48 89 d6             	mov    %rdx,%rsi
    c8c1:	48 89 c7             	mov    %rax,%rdi
    c8c4:	e8 17 e2 ff ff       	call   aae0 <MPI_Type_set_name@plt>
}
    c8c9:	90                   	nop
    c8ca:	c9                   	leave  
    c8cb:	c3                   	ret    

000000000000c8cc <_ZN3MPI14Is_initializedEv>:
  (void)MPI_Finalize();
}

inline bool
MPI::Is_initialized()
{
    c8cc:	f3 0f 1e fa          	endbr64 
    c8d0:	55                   	push   %rbp
    c8d1:	48 89 e5             	mov    %rsp,%rbp
    c8d4:	48 83 ec 10          	sub    $0x10,%rsp
    c8d8:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c8df:	00 00 
    c8e1:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c8e5:	31 c0                	xor    %eax,%eax
  int t;
  (void)MPI_Initialized(&t);
    c8e7:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
    c8eb:	48 89 c7             	mov    %rax,%rdi
    c8ee:	e8 8d d7 ff ff       	call   a080 <MPI_Initialized@plt>
  return OPAL_INT_TO_BOOL(t);
    c8f3:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c8f6:	85 c0                	test   %eax,%eax
    c8f8:	0f 95 c0             	setne  %al
}
    c8fb:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    c8ff:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    c906:	00 00 
    c908:	74 05                	je     c90f <_ZN3MPI14Is_initializedEv+0x43>
    c90a:	e8 d1 db ff ff       	call   a4e0 <__stack_chk_fail@plt>
    c90f:	c9                   	leave  
    c910:	c3                   	ret    
    c911:	90                   	nop

000000000000c912 <_ZN3MPI7Request4WaitERNS_6StatusE>:
// Point-to-Point Communication
//

inline void
MPI::Request::Wait(MPI::Status &status)
{
    c912:	f3 0f 1e fa          	endbr64 
    c916:	55                   	push   %rbp
    c917:	48 89 e5             	mov    %rsp,%rbp
    c91a:	48 83 ec 10          	sub    $0x10,%rsp
    c91e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c922:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
  (void)MPI_Wait(&mpi_request, &status.mpi_status);
    c926:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    c92a:	48 8d 50 08          	lea    0x8(%rax),%rdx
    c92e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c932:	48 83 c0 08          	add    $0x8,%rax
    c936:	48 89 d6             	mov    %rdx,%rsi
    c939:	48 89 c7             	mov    %rax,%rdi
    c93c:	e8 ff dc ff ff       	call   a640 <MPI_Wait@plt>
}
    c941:	90                   	nop
    c942:	c9                   	leave  
    c943:	c3                   	ret    

000000000000c944 <_ZN3MPI7Request4WaitEv>:

inline void
MPI::Request::Wait()
{
    c944:	f3 0f 1e fa          	endbr64 
    c948:	55                   	push   %rbp
    c949:	48 89 e5             	mov    %rsp,%rbp
    c94c:	48 83 ec 10          	sub    $0x10,%rsp
    c950:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
  (void)MPI_Wait(&mpi_request, MPI_STATUS_IGNORE);
    c954:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c958:	48 83 c0 08          	add    $0x8,%rax
    c95c:	be 00 00 00 00       	mov    $0x0,%esi
    c961:	48 89 c7             	mov    %rax,%rdi
    c964:	e8 d7 dc ff ff       	call   a640 <MPI_Wait@plt>
}
    c969:	90                   	nop
    c96a:	c9                   	leave  
    c96b:	c3                   	ret    

000000000000c96c <_ZN3MPI7Request4FreeEv>:

inline void
MPI::Request::Free()
{
    c96c:	f3 0f 1e fa          	endbr64 
    c970:	55                   	push   %rbp
    c971:	48 89 e5             	mov    %rsp,%rbp
    c974:	48 83 ec 10          	sub    $0x10,%rsp
    c978:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
  (void)MPI_Request_free(&mpi_request);
    c97c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c980:	48 83 c0 08          	add    $0x8,%rax
    c984:	48 89 c7             	mov    %rax,%rdi
    c987:	e8 f4 d8 ff ff       	call   a280 <MPI_Request_free@plt>
}
    c98c:	90                   	nop
    c98d:	c9                   	leave  
    c98e:	c3                   	ret    
    c98f:	90                   	nop

000000000000c990 <_ZN3MPI7Request4TestERNS_6StatusE>:

inline bool
MPI::Request::Test(MPI::Status &status)
{
    c990:	f3 0f 1e fa          	endbr64 
    c994:	55                   	push   %rbp
    c995:	48 89 e5             	mov    %rsp,%rbp
    c998:	48 83 ec 20          	sub    $0x20,%rsp
    c99c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c9a0:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c9a4:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c9ab:	00 00 
    c9ad:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c9b1:	31 c0                	xor    %eax,%eax
  int t;
  (void)MPI_Test(&mpi_request, &t, &status.mpi_status);
    c9b3:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c9b7:	48 8d 50 08          	lea    0x8(%rax),%rdx
    c9bb:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c9bf:	48 8d 48 08          	lea    0x8(%rax),%rcx
    c9c3:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
    c9c7:	48 89 c6             	mov    %rax,%rsi
    c9ca:	48 89 cf             	mov    %rcx,%rdi
    c9cd:	e8 0e df ff ff       	call   a8e0 <MPI_Test@plt>
  return OPAL_INT_TO_BOOL(t);
    c9d2:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c9d5:	85 c0                	test   %eax,%eax
    c9d7:	0f 95 c0             	setne  %al
}
    c9da:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    c9de:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    c9e5:	00 00 
    c9e7:	74 05                	je     c9ee <_ZN3MPI7Request4TestERNS_6StatusE+0x5e>
    c9e9:	e8 f2 da ff ff       	call   a4e0 <__stack_chk_fail@plt>
    c9ee:	c9                   	leave  
    c9ef:	c3                   	ret    

000000000000c9f0 <_ZN3MPI7Request4TestEv>:

inline bool
MPI::Request::Test()
{
    c9f0:	f3 0f 1e fa          	endbr64 
    c9f4:	55                   	push   %rbp
    c9f5:	48 89 e5             	mov    %rsp,%rbp
    c9f8:	48 83 ec 20          	sub    $0x20,%rsp
    c9fc:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    ca00:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ca07:	00 00 
    ca09:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    ca0d:	31 c0                	xor    %eax,%eax
  int t;
  (void)MPI_Test(&mpi_request, &t, MPI_STATUS_IGNORE);
    ca0f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ca13:	48 8d 48 08          	lea    0x8(%rax),%rcx
    ca17:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
    ca1b:	ba 00 00 00 00       	mov    $0x0,%edx
    ca20:	48 89 c6             	mov    %rax,%rsi
    ca23:	48 89 cf             	mov    %rcx,%rdi
    ca26:	e8 b5 de ff ff       	call   a8e0 <MPI_Test@plt>
  return OPAL_INT_TO_BOOL(t);
    ca2b:	8b 45 f4             	mov    -0xc(%rbp),%eax
    ca2e:	85 c0                	test   %eax,%eax
    ca30:	0f 95 c0             	setne  %al
}
    ca33:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    ca37:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    ca3e:	00 00 
    ca40:	74 05                	je     ca47 <_ZN3MPI7Request4TestEv+0x57>
    ca42:	e8 99 da ff ff       	call   a4e0 <__stack_chk_fail@plt>
    ca47:	c9                   	leave  
    ca48:	c3                   	ret    
    ca49:	90                   	nop

000000000000ca4a <_ZNK3MPI7Request6CancelEv>:
  return outcount;
}

inline void
MPI::Request::Cancel(void) const
{
    ca4a:	f3 0f 1e fa          	endbr64 
    ca4e:	55                   	push   %rbp
    ca4f:	48 89 e5             	mov    %rsp,%rbp
    ca52:	48 83 ec 10          	sub    $0x10,%rsp
    ca56:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
  (void)MPI_Cancel(const_cast<MPI_Request *>(&mpi_request));
    ca5a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ca5e:	48 83 c0 08          	add    $0x8,%rax
    ca62:	48 89 c7             	mov    %rax,%rdi
    ca65:	e8 06 d6 ff ff       	call   a070 <MPI_Cancel@plt>
}
    ca6a:	90                   	nop
    ca6b:	c9                   	leave  
    ca6c:	c3                   	ret    
    ca6d:	90                   	nop

000000000000ca6e <_ZN3MPI8Prequest5StartEv>:

inline void
MPI::Prequest::Start()
{
    ca6e:	f3 0f 1e fa          	endbr64 
    ca72:	55                   	push   %rbp
    ca73:	48 89 e5             	mov    %rsp,%rbp
    ca76:	48 83 ec 10          	sub    $0x10,%rsp
    ca7a:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
  (void)MPI_Start(&mpi_request);
    ca7e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ca82:	48 83 c0 08          	add    $0x8,%rax
    ca86:	48 89 c7             	mov    %rax,%rdi
    ca89:	e8 92 de ff ff       	call   a920 <MPI_Start@plt>
}
    ca8e:	90                   	nop
    ca8f:	c9                   	leave  
    ca90:	c3                   	ret    
    ca91:	90                   	nop

000000000000ca92 <_ZNK3MPI7Request10Get_statusERNS_6StatusE>:
  }
  delete [] mpi_requests;
}

inline bool MPI::Request::Get_status(MPI::Status& status) const
{
    ca92:	f3 0f 1e fa          	endbr64 
    ca96:	55                   	push   %rbp
    ca97:	48 89 e5             	mov    %rsp,%rbp
    ca9a:	48 83 ec 40          	sub    $0x40,%rsp
    ca9e:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
    caa2:	48 89 75 c0          	mov    %rsi,-0x40(%rbp)
    caa6:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    caad:	00 00 
    caaf:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    cab3:	31 c0                	xor    %eax,%eax
    int flag = 0;
    cab5:	c7 45 dc 00 00 00 00 	movl   $0x0,-0x24(%rbp)
    MPI_Status c_status;

    // Call the underlying MPI function rather than simply returning
    // status.mpi_status because we may have to invoke the generalized
    // request query function
    (void)MPI_Request_get_status(mpi_request, &flag, &c_status);
    cabc:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    cac0:	48 8b 40 08          	mov    0x8(%rax),%rax
    cac4:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    cac8:	48 8d 4d dc          	lea    -0x24(%rbp),%rcx
    cacc:	48 89 ce             	mov    %rcx,%rsi
    cacf:	48 89 c7             	mov    %rax,%rdi
    cad2:	e8 69 d4 ff ff       	call   9f40 <MPI_Request_get_status@plt>
    if (flag) {
    cad7:	8b 45 dc             	mov    -0x24(%rbp),%eax
    cada:	85 c0                	test   %eax,%eax
    cadc:	74 13                	je     caf1 <_ZNK3MPI7Request10Get_statusERNS_6StatusE+0x5f>
        status = c_status;
    cade:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    cae2:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    cae6:	48 89 d6             	mov    %rdx,%rsi
    cae9:	48 89 c7             	mov    %rax,%rdi
    caec:	e8 f1 ed ff ff       	call   b8e2 <_ZN3MPI6StatusaSERK20ompi_status_public_t>
    }
    return OPAL_INT_TO_BOOL(flag);
    caf1:	8b 45 dc             	mov    -0x24(%rbp),%eax
    caf4:	85 c0                	test   %eax,%eax
    caf6:	0f 95 c0             	setne  %al
}
    caf9:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    cafd:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    cb04:	00 00 
    cb06:	74 05                	je     cb0d <_ZNK3MPI7Request10Get_statusERNS_6StatusE+0x7b>
    cb08:	e8 d3 d9 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    cb0d:	c9                   	leave  
    cb0e:	c3                   	ret    
    cb0f:	90                   	nop

000000000000cb10 <_ZNK3MPI7Request10Get_statusEv>:

inline bool MPI::Request::Get_status() const
{
    cb10:	f3 0f 1e fa          	endbr64 
    cb14:	55                   	push   %rbp
    cb15:	48 89 e5             	mov    %rsp,%rbp
    cb18:	48 83 ec 20          	sub    $0x20,%rsp
    cb1c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cb20:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    cb27:	00 00 
    cb29:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    cb2d:	31 c0                	xor    %eax,%eax
    int flag;

    // Call the underlying MPI function rather than simply returning
    // status.mpi_status because we may have to invoke the generalized
    // request query function
    (void)MPI_Request_get_status(mpi_request, &flag, MPI_STATUS_IGNORE);
    cb2f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cb33:	48 8b 40 08          	mov    0x8(%rax),%rax
    cb37:	48 8d 4d f4          	lea    -0xc(%rbp),%rcx
    cb3b:	ba 00 00 00 00       	mov    $0x0,%edx
    cb40:	48 89 ce             	mov    %rcx,%rsi
    cb43:	48 89 c7             	mov    %rax,%rdi
    cb46:	e8 f5 d3 ff ff       	call   9f40 <MPI_Request_get_status@plt>
    return OPAL_INT_TO_BOOL(flag);
    cb4b:	8b 45 f4             	mov    -0xc(%rbp),%eax
    cb4e:	85 c0                	test   %eax,%eax
    cb50:	0f 95 c0             	setne  %al
}
    cb53:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    cb57:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    cb5e:	00 00 
    cb60:	74 05                	je     cb67 <_ZNK3MPI7Request10Get_statusEv+0x57>
    cb62:	e8 79 d9 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    cb67:	c9                   	leave  
    cb68:	c3                   	ret    
    cb69:	90                   	nop

000000000000cb6a <_ZN3MPI8Grequest8CompleteEv>:
    return(grequest);
}

inline void
MPI::Grequest::Complete()
{
    cb6a:	f3 0f 1e fa          	endbr64 
    cb6e:	55                   	push   %rbp
    cb6f:	48 89 e5             	mov    %rsp,%rbp
    cb72:	48 83 ec 10          	sub    $0x10,%rsp
    cb76:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    (void) MPI_Grequest_complete(mpi_request);
    cb7a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    cb7e:	48 8b 40 08          	mov    0x8(%rax),%rax
    cb82:	48 89 c7             	mov    %rax,%rdi
    cb85:	e8 16 d7 ff ff       	call   a2a0 <MPI_Grequest_complete@plt>
}
    cb8a:	90                   	nop
    cb8b:	c9                   	leave  
    cb8c:	c3                   	ret    
    cb8d:	90                   	nop

000000000000cb8e <_ZNK3MPI4Comm4SendEPKviRKNS_8DatatypeEii>:
//

inline void
MPI::Comm::Send(const void *buf, int count,
		const MPI::Datatype & datatype, int dest, int tag) const
{
    cb8e:	f3 0f 1e fa          	endbr64 
    cb92:	55                   	push   %rbp
    cb93:	48 89 e5             	mov    %rsp,%rbp
    cb96:	53                   	push   %rbx
    cb97:	48 83 ec 38          	sub    $0x38,%rsp
    cb9b:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cb9f:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    cba3:	89 55 dc             	mov    %edx,-0x24(%rbp)
    cba6:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    cbaa:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    cbae:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
  (void)MPI_Send(const_cast<void *>(buf), count, datatype, dest, tag, mpi_comm);
    cbb2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cbb6:	48 8b 58 08          	mov    0x8(%rax),%rbx
    cbba:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cbbe:	48 89 c7             	mov    %rax,%rdi
    cbc1:	e8 b8 ec ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cbc6:	48 89 c7             	mov    %rax,%rdi
    cbc9:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    cbcc:	8b 55 d8             	mov    -0x28(%rbp),%edx
    cbcf:	8b 75 dc             	mov    -0x24(%rbp),%esi
    cbd2:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    cbd6:	49 89 d9             	mov    %rbx,%r9
    cbd9:	41 89 c8             	mov    %ecx,%r8d
    cbdc:	89 d1                	mov    %edx,%ecx
    cbde:	48 89 fa             	mov    %rdi,%rdx
    cbe1:	48 89 c7             	mov    %rax,%rdi
    cbe4:	e8 b7 d3 ff ff       	call   9fa0 <MPI_Send@plt>
}
    cbe9:	90                   	nop
    cbea:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    cbee:	c9                   	leave  
    cbef:	c3                   	ret    

000000000000cbf0 <_ZNK3MPI4Comm4RecvEPviRKNS_8DatatypeEiiRNS_6StatusE>:

inline void
MPI::Comm::Recv(void *buf, int count, const MPI::Datatype & datatype,
		int source, int tag, MPI::Status & status) const
{
    cbf0:	f3 0f 1e fa          	endbr64 
    cbf4:	55                   	push   %rbp
    cbf5:	48 89 e5             	mov    %rsp,%rbp
    cbf8:	41 54                	push   %r12
    cbfa:	53                   	push   %rbx
    cbfb:	48 83 ec 30          	sub    $0x30,%rsp
    cbff:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cc03:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    cc07:	89 55 dc             	mov    %edx,-0x24(%rbp)
    cc0a:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    cc0e:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    cc12:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
  (void)MPI_Recv(buf, count, datatype, source, tag, mpi_comm, &status.mpi_status);
    cc16:	48 8b 45 10          	mov    0x10(%rbp),%rax
    cc1a:	48 8d 58 08          	lea    0x8(%rax),%rbx
    cc1e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cc22:	4c 8b 60 08          	mov    0x8(%rax),%r12
    cc26:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cc2a:	48 89 c7             	mov    %rax,%rdi
    cc2d:	e8 4c ec ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cc32:	48 89 c7             	mov    %rax,%rdi
    cc35:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    cc38:	8b 55 d8             	mov    -0x28(%rbp),%edx
    cc3b:	8b 75 dc             	mov    -0x24(%rbp),%esi
    cc3e:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    cc42:	48 83 ec 08          	sub    $0x8,%rsp
    cc46:	53                   	push   %rbx
    cc47:	4d 89 e1             	mov    %r12,%r9
    cc4a:	41 89 c8             	mov    %ecx,%r8d
    cc4d:	89 d1                	mov    %edx,%ecx
    cc4f:	48 89 fa             	mov    %rdi,%rdx
    cc52:	48 89 c7             	mov    %rax,%rdi
    cc55:	e8 86 d6 ff ff       	call   a2e0 <MPI_Recv@plt>
    cc5a:	48 83 c4 10          	add    $0x10,%rsp
}
    cc5e:	90                   	nop
    cc5f:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    cc63:	5b                   	pop    %rbx
    cc64:	41 5c                	pop    %r12
    cc66:	5d                   	pop    %rbp
    cc67:	c3                   	ret    

000000000000cc68 <_ZNK3MPI4Comm4RecvEPviRKNS_8DatatypeEii>:

inline void
MPI::Comm::Recv(void *buf, int count, const MPI::Datatype & datatype,
				    int source, int tag) const
{
    cc68:	f3 0f 1e fa          	endbr64 
    cc6c:	55                   	push   %rbp
    cc6d:	48 89 e5             	mov    %rsp,%rbp
    cc70:	53                   	push   %rbx
    cc71:	48 83 ec 38          	sub    $0x38,%rsp
    cc75:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cc79:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    cc7d:	89 55 dc             	mov    %edx,-0x24(%rbp)
    cc80:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    cc84:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    cc88:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
  (void)MPI_Recv(buf, count, datatype, source,
		 tag, mpi_comm, MPI_STATUS_IGNORE);
    cc8c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cc90:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Recv(buf, count, datatype, source,
    cc94:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cc98:	48 89 c7             	mov    %rax,%rdi
    cc9b:	e8 de eb ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cca0:	48 89 c7             	mov    %rax,%rdi
    cca3:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    cca6:	8b 55 d8             	mov    -0x28(%rbp),%edx
    cca9:	8b 75 dc             	mov    -0x24(%rbp),%esi
    ccac:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    ccb0:	48 83 ec 08          	sub    $0x8,%rsp
    ccb4:	6a 00                	push   $0x0
    ccb6:	49 89 d9             	mov    %rbx,%r9
    ccb9:	41 89 c8             	mov    %ecx,%r8d
    ccbc:	89 d1                	mov    %edx,%ecx
    ccbe:	48 89 fa             	mov    %rdi,%rdx
    ccc1:	48 89 c7             	mov    %rax,%rdi
    ccc4:	e8 17 d6 ff ff       	call   a2e0 <MPI_Recv@plt>
    ccc9:	48 83 c4 10          	add    $0x10,%rsp
}
    cccd:	90                   	nop
    ccce:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    ccd2:	c9                   	leave  
    ccd3:	c3                   	ret    

000000000000ccd4 <_ZNK3MPI4Comm5BsendEPKviRKNS_8DatatypeEii>:

inline void
MPI::Comm::Bsend(const void *buf, int count,
		 const MPI::Datatype & datatype, int dest, int tag) const
{
    ccd4:	f3 0f 1e fa          	endbr64 
    ccd8:	55                   	push   %rbp
    ccd9:	48 89 e5             	mov    %rsp,%rbp
    ccdc:	53                   	push   %rbx
    ccdd:	48 83 ec 38          	sub    $0x38,%rsp
    cce1:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cce5:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    cce9:	89 55 dc             	mov    %edx,-0x24(%rbp)
    ccec:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    ccf0:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    ccf4:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
  (void)MPI_Bsend(const_cast<void *>(buf), count, datatype,
		  dest, tag, mpi_comm);
    ccf8:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ccfc:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Bsend(const_cast<void *>(buf), count, datatype,
    cd00:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cd04:	48 89 c7             	mov    %rax,%rdi
    cd07:	e8 72 eb ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cd0c:	48 89 c7             	mov    %rax,%rdi
    cd0f:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    cd12:	8b 55 d8             	mov    -0x28(%rbp),%edx
    cd15:	8b 75 dc             	mov    -0x24(%rbp),%esi
    cd18:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    cd1c:	49 89 d9             	mov    %rbx,%r9
    cd1f:	41 89 c8             	mov    %ecx,%r8d
    cd22:	89 d1                	mov    %edx,%ecx
    cd24:	48 89 fa             	mov    %rdi,%rdx
    cd27:	48 89 c7             	mov    %rax,%rdi
    cd2a:	e8 91 d7 ff ff       	call   a4c0 <MPI_Bsend@plt>
}
    cd2f:	90                   	nop
    cd30:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    cd34:	c9                   	leave  
    cd35:	c3                   	ret    

000000000000cd36 <_ZNK3MPI4Comm5SsendEPKviRKNS_8DatatypeEii>:

inline void
MPI::Comm::Ssend(const void *buf, int count,
		 const MPI::Datatype & datatype, int dest, int tag) const
{
    cd36:	f3 0f 1e fa          	endbr64 
    cd3a:	55                   	push   %rbp
    cd3b:	48 89 e5             	mov    %rsp,%rbp
    cd3e:	53                   	push   %rbx
    cd3f:	48 83 ec 38          	sub    $0x38,%rsp
    cd43:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cd47:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    cd4b:	89 55 dc             	mov    %edx,-0x24(%rbp)
    cd4e:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    cd52:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    cd56:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
  (void)MPI_Ssend(const_cast<void *>(buf), count,  datatype, dest,
		  tag, mpi_comm);
    cd5a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cd5e:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Ssend(const_cast<void *>(buf), count,  datatype, dest,
    cd62:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cd66:	48 89 c7             	mov    %rax,%rdi
    cd69:	e8 10 eb ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cd6e:	48 89 c7             	mov    %rax,%rdi
    cd71:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    cd74:	8b 55 d8             	mov    -0x28(%rbp),%edx
    cd77:	8b 75 dc             	mov    -0x24(%rbp),%esi
    cd7a:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    cd7e:	49 89 d9             	mov    %rbx,%r9
    cd81:	41 89 c8             	mov    %ecx,%r8d
    cd84:	89 d1                	mov    %edx,%ecx
    cd86:	48 89 fa             	mov    %rdi,%rdx
    cd89:	48 89 c7             	mov    %rax,%rdi
    cd8c:	e8 ff d0 ff ff       	call   9e90 <MPI_Ssend@plt>
}
    cd91:	90                   	nop
    cd92:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    cd96:	c9                   	leave  
    cd97:	c3                   	ret    

000000000000cd98 <_ZNK3MPI4Comm5RsendEPKviRKNS_8DatatypeEii>:

inline void
MPI::Comm::Rsend(const void *buf, int count,
		 const MPI::Datatype & datatype, int dest, int tag) const
{
    cd98:	f3 0f 1e fa          	endbr64 
    cd9c:	55                   	push   %rbp
    cd9d:	48 89 e5             	mov    %rsp,%rbp
    cda0:	53                   	push   %rbx
    cda1:	48 83 ec 38          	sub    $0x38,%rsp
    cda5:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cda9:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    cdad:	89 55 dc             	mov    %edx,-0x24(%rbp)
    cdb0:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    cdb4:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    cdb8:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
  (void)MPI_Rsend(const_cast<void *>(buf), count, datatype,
		  dest, tag, mpi_comm);
    cdbc:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cdc0:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Rsend(const_cast<void *>(buf), count, datatype,
    cdc4:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cdc8:	48 89 c7             	mov    %rax,%rdi
    cdcb:	e8 ae ea ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cdd0:	48 89 c7             	mov    %rax,%rdi
    cdd3:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    cdd6:	8b 55 d8             	mov    -0x28(%rbp),%edx
    cdd9:	8b 75 dc             	mov    -0x24(%rbp),%esi
    cddc:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    cde0:	49 89 d9             	mov    %rbx,%r9
    cde3:	41 89 c8             	mov    %ecx,%r8d
    cde6:	89 d1                	mov    %edx,%ecx
    cde8:	48 89 fa             	mov    %rdi,%rdx
    cdeb:	48 89 c7             	mov    %rax,%rdi
    cdee:	e8 dd d4 ff ff       	call   a2d0 <MPI_Rsend@plt>
}
    cdf3:	90                   	nop
    cdf4:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    cdf8:	c9                   	leave  
    cdf9:	c3                   	ret    

000000000000cdfa <_ZNK3MPI4Comm5IsendEPKviRKNS_8DatatypeEii>:

inline MPI::Request
MPI::Comm::Isend(const void *buf, int count,
		 const MPI::Datatype & datatype, int dest, int tag) const
{
    cdfa:	f3 0f 1e fa          	endbr64 
    cdfe:	55                   	push   %rbp
    cdff:	48 89 e5             	mov    %rsp,%rbp
    ce02:	53                   	push   %rbx
    ce03:	48 83 ec 48          	sub    $0x48,%rsp
    ce07:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    ce0b:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    ce0f:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    ce13:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    ce16:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    ce1a:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    ce1e:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ce25:	00 00 
    ce27:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    ce2b:	31 c0                	xor    %eax,%eax
  MPI_Request request;
  (void)MPI_Isend(const_cast<void *>(buf), count, datatype,
		  dest, tag, mpi_comm, &request);
    ce2d:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ce31:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Isend(const_cast<void *>(buf), count, datatype,
    ce35:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    ce39:	48 89 c7             	mov    %rax,%rdi
    ce3c:	e8 3d ea ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    ce41:	48 89 c7             	mov    %rax,%rdi
    ce44:	8b 55 c0             	mov    -0x40(%rbp),%edx
    ce47:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    ce4a:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    ce4e:	48 83 ec 08          	sub    $0x8,%rsp
    ce52:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    ce56:	51                   	push   %rcx
    ce57:	49 89 d9             	mov    %rbx,%r9
    ce5a:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    ce5e:	89 d1                	mov    %edx,%ecx
    ce60:	48 89 fa             	mov    %rdi,%rdx
    ce63:	48 89 c7             	mov    %rax,%rdi
    ce66:	e8 a5 d1 ff ff       	call   a010 <MPI_Isend@plt>
    ce6b:	48 83 c4 10          	add    $0x10,%rsp
  return request;
    ce6f:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    ce73:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ce77:	48 89 d6             	mov    %rdx,%rsi
    ce7a:	48 89 c7             	mov    %rax,%rdi
    ce7d:	e8 e4 ea ff ff       	call   b966 <_ZN3MPI7RequestC1EP14ompi_request_t>
}
    ce82:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ce86:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    ce8d:	00 00 
    ce8f:	74 05                	je     ce96 <_ZNK3MPI4Comm5IsendEPKviRKNS_8DatatypeEii+0x9c>
    ce91:	e8 4a d6 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    ce96:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ce9a:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    ce9e:	c9                   	leave  
    ce9f:	c3                   	ret    

000000000000cea0 <_ZNK3MPI4Comm6IbsendEPKviRKNS_8DatatypeEii>:

inline MPI::Request
MPI::Comm::Ibsend(const void *buf, int count,
		  const MPI::Datatype & datatype, int dest, int tag) const
{
    cea0:	f3 0f 1e fa          	endbr64 
    cea4:	55                   	push   %rbp
    cea5:	48 89 e5             	mov    %rsp,%rbp
    cea8:	53                   	push   %rbx
    cea9:	48 83 ec 48          	sub    $0x48,%rsp
    cead:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    ceb1:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    ceb5:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    ceb9:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    cebc:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    cec0:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    cec4:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    cecb:	00 00 
    cecd:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    ced1:	31 c0                	xor    %eax,%eax
  MPI_Request request;
  (void)MPI_Ibsend(const_cast<void *>(buf), count, datatype,
		   dest, tag, mpi_comm, &request);
    ced3:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ced7:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Ibsend(const_cast<void *>(buf), count, datatype,
    cedb:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    cedf:	48 89 c7             	mov    %rax,%rdi
    cee2:	e8 97 e9 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cee7:	48 89 c7             	mov    %rax,%rdi
    ceea:	8b 55 c0             	mov    -0x40(%rbp),%edx
    ceed:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    cef0:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    cef4:	48 83 ec 08          	sub    $0x8,%rsp
    cef8:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    cefc:	51                   	push   %rcx
    cefd:	49 89 d9             	mov    %rbx,%r9
    cf00:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    cf04:	89 d1                	mov    %edx,%ecx
    cf06:	48 89 fa             	mov    %rdi,%rdx
    cf09:	48 89 c7             	mov    %rax,%rdi
    cf0c:	e8 ef d8 ff ff       	call   a800 <MPI_Ibsend@plt>
    cf11:	48 83 c4 10          	add    $0x10,%rsp
  return request;
    cf15:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    cf19:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    cf1d:	48 89 d6             	mov    %rdx,%rsi
    cf20:	48 89 c7             	mov    %rax,%rdi
    cf23:	e8 3e ea ff ff       	call   b966 <_ZN3MPI7RequestC1EP14ompi_request_t>
}
    cf28:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cf2c:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    cf33:	00 00 
    cf35:	74 05                	je     cf3c <_ZNK3MPI4Comm6IbsendEPKviRKNS_8DatatypeEii+0x9c>
    cf37:	e8 a4 d5 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    cf3c:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    cf40:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    cf44:	c9                   	leave  
    cf45:	c3                   	ret    

000000000000cf46 <_ZNK3MPI4Comm6IssendEPKviRKNS_8DatatypeEii>:

inline MPI::Request
MPI::Comm::Issend(const void *buf, int count,
		  const MPI::Datatype & datatype, int dest, int tag) const
{
    cf46:	f3 0f 1e fa          	endbr64 
    cf4a:	55                   	push   %rbp
    cf4b:	48 89 e5             	mov    %rsp,%rbp
    cf4e:	53                   	push   %rbx
    cf4f:	48 83 ec 48          	sub    $0x48,%rsp
    cf53:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    cf57:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    cf5b:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    cf5f:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    cf62:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    cf66:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    cf6a:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    cf71:	00 00 
    cf73:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    cf77:	31 c0                	xor    %eax,%eax
  MPI_Request request;
  (void)MPI_Issend(const_cast<void *>(buf), count, datatype,
		   dest, tag, mpi_comm, &request);
    cf79:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cf7d:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Issend(const_cast<void *>(buf), count, datatype,
    cf81:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    cf85:	48 89 c7             	mov    %rax,%rdi
    cf88:	e8 f1 e8 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cf8d:	48 89 c7             	mov    %rax,%rdi
    cf90:	8b 55 c0             	mov    -0x40(%rbp),%edx
    cf93:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    cf96:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    cf9a:	48 83 ec 08          	sub    $0x8,%rsp
    cf9e:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    cfa2:	51                   	push   %rcx
    cfa3:	49 89 d9             	mov    %rbx,%r9
    cfa6:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    cfaa:	89 d1                	mov    %edx,%ecx
    cfac:	48 89 fa             	mov    %rdi,%rdx
    cfaf:	48 89 c7             	mov    %rax,%rdi
    cfb2:	e8 e9 d0 ff ff       	call   a0a0 <MPI_Issend@plt>
    cfb7:	48 83 c4 10          	add    $0x10,%rsp
  return request;
    cfbb:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    cfbf:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    cfc3:	48 89 d6             	mov    %rdx,%rsi
    cfc6:	48 89 c7             	mov    %rax,%rdi
    cfc9:	e8 98 e9 ff ff       	call   b966 <_ZN3MPI7RequestC1EP14ompi_request_t>
}
    cfce:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cfd2:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    cfd9:	00 00 
    cfdb:	74 05                	je     cfe2 <_ZNK3MPI4Comm6IssendEPKviRKNS_8DatatypeEii+0x9c>
    cfdd:	e8 fe d4 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    cfe2:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    cfe6:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    cfea:	c9                   	leave  
    cfeb:	c3                   	ret    

000000000000cfec <_ZNK3MPI4Comm6IrsendEPKviRKNS_8DatatypeEii>:

inline MPI::Request
MPI::Comm::Irsend(const void *buf, int count,
		  const MPI::Datatype & datatype, int dest, int tag) const
{
    cfec:	f3 0f 1e fa          	endbr64 
    cff0:	55                   	push   %rbp
    cff1:	48 89 e5             	mov    %rsp,%rbp
    cff4:	53                   	push   %rbx
    cff5:	48 83 ec 48          	sub    $0x48,%rsp
    cff9:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    cffd:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d001:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d005:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d008:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d00c:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d010:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d017:	00 00 
    d019:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d01d:	31 c0                	xor    %eax,%eax
  MPI_Request request;
  (void)MPI_Irsend(const_cast<void *>(buf), count, datatype,
		   dest, tag, mpi_comm, &request);
    d01f:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d023:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Irsend(const_cast<void *>(buf), count, datatype,
    d027:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d02b:	48 89 c7             	mov    %rax,%rdi
    d02e:	e8 4b e8 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d033:	48 89 c7             	mov    %rax,%rdi
    d036:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d039:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d03c:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d040:	48 83 ec 08          	sub    $0x8,%rsp
    d044:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d048:	51                   	push   %rcx
    d049:	49 89 d9             	mov    %rbx,%r9
    d04c:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d050:	89 d1                	mov    %edx,%ecx
    d052:	48 89 fa             	mov    %rdi,%rdx
    d055:	48 89 c7             	mov    %rax,%rdi
    d058:	e8 03 ce ff ff       	call   9e60 <MPI_Irsend@plt>
    d05d:	48 83 c4 10          	add    $0x10,%rsp
  return request;
    d061:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    d065:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d069:	48 89 d6             	mov    %rdx,%rsi
    d06c:	48 89 c7             	mov    %rax,%rdi
    d06f:	e8 f2 e8 ff ff       	call   b966 <_ZN3MPI7RequestC1EP14ompi_request_t>
}
    d074:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d078:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d07f:	00 00 
    d081:	74 05                	je     d088 <_ZNK3MPI4Comm6IrsendEPKviRKNS_8DatatypeEii+0x9c>
    d083:	e8 58 d4 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    d088:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d08c:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d090:	c9                   	leave  
    d091:	c3                   	ret    

000000000000d092 <_ZNK3MPI4Comm5IrecvEPviRKNS_8DatatypeEii>:

inline MPI::Request
MPI::Comm::Irecv(void *buf, int count,
		 const MPI::Datatype & datatype, int source, int tag) const
{
    d092:	f3 0f 1e fa          	endbr64 
    d096:	55                   	push   %rbp
    d097:	48 89 e5             	mov    %rsp,%rbp
    d09a:	53                   	push   %rbx
    d09b:	48 83 ec 48          	sub    $0x48,%rsp
    d09f:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d0a3:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d0a7:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d0ab:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d0ae:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d0b2:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d0b6:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d0bd:	00 00 
    d0bf:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d0c3:	31 c0                	xor    %eax,%eax
  MPI_Request request;
  (void)MPI_Irecv(buf, count, datatype, source,
		  tag, mpi_comm, &request);
    d0c5:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d0c9:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Irecv(buf, count, datatype, source,
    d0cd:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d0d1:	48 89 c7             	mov    %rax,%rdi
    d0d4:	e8 a5 e7 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d0d9:	48 89 c7             	mov    %rax,%rdi
    d0dc:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d0df:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d0e2:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d0e6:	48 83 ec 08          	sub    $0x8,%rsp
    d0ea:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d0ee:	51                   	push   %rcx
    d0ef:	49 89 d9             	mov    %rbx,%r9
    d0f2:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d0f6:	89 d1                	mov    %edx,%ecx
    d0f8:	48 89 fa             	mov    %rdi,%rdx
    d0fb:	48 89 c7             	mov    %rax,%rdi
    d0fe:	e8 7d d5 ff ff       	call   a680 <MPI_Irecv@plt>
    d103:	48 83 c4 10          	add    $0x10,%rsp
  return request;
    d107:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    d10b:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d10f:	48 89 d6             	mov    %rdx,%rsi
    d112:	48 89 c7             	mov    %rax,%rdi
    d115:	e8 4c e8 ff ff       	call   b966 <_ZN3MPI7RequestC1EP14ompi_request_t>
}
    d11a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d11e:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d125:	00 00 
    d127:	74 05                	je     d12e <_ZNK3MPI4Comm5IrecvEPviRKNS_8DatatypeEii+0x9c>
    d129:	e8 b2 d3 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    d12e:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d132:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d136:	c9                   	leave  
    d137:	c3                   	ret    

000000000000d138 <_ZNK3MPI4Comm6IprobeEiiRNS_6StatusE>:


inline bool
MPI::Comm::Iprobe(int source, int tag, MPI::Status & status) const
{
    d138:	f3 0f 1e fa          	endbr64 
    d13c:	55                   	push   %rbp
    d13d:	48 89 e5             	mov    %rsp,%rbp
    d140:	48 83 ec 30          	sub    $0x30,%rsp
    d144:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d148:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    d14b:	89 55 e0             	mov    %edx,-0x20(%rbp)
    d14e:	48 89 4d d8          	mov    %rcx,-0x28(%rbp)
    d152:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d159:	00 00 
    d15b:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    d15f:	31 c0                	xor    %eax,%eax
  int t;
  (void)MPI_Iprobe(source, tag, mpi_comm, &t, &status.mpi_status);
    d161:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d165:	48 8d 78 08          	lea    0x8(%rax),%rdi
    d169:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d16d:	48 8b 50 08          	mov    0x8(%rax),%rdx
    d171:	48 8d 4d f4          	lea    -0xc(%rbp),%rcx
    d175:	8b 75 e0             	mov    -0x20(%rbp),%esi
    d178:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    d17b:	49 89 f8             	mov    %rdi,%r8
    d17e:	89 c7                	mov    %eax,%edi
    d180:	e8 ab ce ff ff       	call   a030 <MPI_Iprobe@plt>
  return OPAL_INT_TO_BOOL(t);
    d185:	8b 45 f4             	mov    -0xc(%rbp),%eax
    d188:	85 c0                	test   %eax,%eax
    d18a:	0f 95 c0             	setne  %al
}
    d18d:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    d191:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    d198:	00 00 
    d19a:	74 05                	je     d1a1 <_ZNK3MPI4Comm6IprobeEiiRNS_6StatusE+0x69>
    d19c:	e8 3f d3 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    d1a1:	c9                   	leave  
    d1a2:	c3                   	ret    
    d1a3:	90                   	nop

000000000000d1a4 <_ZNK3MPI4Comm6IprobeEii>:

inline bool
MPI::Comm::Iprobe(int source, int tag) const
{
    d1a4:	f3 0f 1e fa          	endbr64 
    d1a8:	55                   	push   %rbp
    d1a9:	48 89 e5             	mov    %rsp,%rbp
    d1ac:	48 83 ec 20          	sub    $0x20,%rsp
    d1b0:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d1b4:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    d1b7:	89 55 e0             	mov    %edx,-0x20(%rbp)
    d1ba:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d1c1:	00 00 
    d1c3:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    d1c7:	31 c0                	xor    %eax,%eax
  int t;
  (void)MPI_Iprobe(source, tag, mpi_comm, &t, MPI_STATUS_IGNORE);
    d1c9:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d1cd:	48 8b 50 08          	mov    0x8(%rax),%rdx
    d1d1:	48 8d 4d f4          	lea    -0xc(%rbp),%rcx
    d1d5:	8b 75 e0             	mov    -0x20(%rbp),%esi
    d1d8:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    d1db:	41 b8 00 00 00 00    	mov    $0x0,%r8d
    d1e1:	89 c7                	mov    %eax,%edi
    d1e3:	e8 48 ce ff ff       	call   a030 <MPI_Iprobe@plt>
  return OPAL_INT_TO_BOOL(t);
    d1e8:	8b 45 f4             	mov    -0xc(%rbp),%eax
    d1eb:	85 c0                	test   %eax,%eax
    d1ed:	0f 95 c0             	setne  %al
}
    d1f0:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    d1f4:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    d1fb:	00 00 
    d1fd:	74 05                	je     d204 <_ZNK3MPI4Comm6IprobeEii+0x60>
    d1ff:	e8 dc d2 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    d204:	c9                   	leave  
    d205:	c3                   	ret    

000000000000d206 <_ZNK3MPI4Comm5ProbeEiiRNS_6StatusE>:

inline void
MPI::Comm::Probe(int source, int tag, MPI::Status & status) const
{
    d206:	f3 0f 1e fa          	endbr64 
    d20a:	55                   	push   %rbp
    d20b:	48 89 e5             	mov    %rsp,%rbp
    d20e:	48 83 ec 20          	sub    $0x20,%rsp
    d212:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    d216:	89 75 f4             	mov    %esi,-0xc(%rbp)
    d219:	89 55 f0             	mov    %edx,-0x10(%rbp)
    d21c:	48 89 4d e8          	mov    %rcx,-0x18(%rbp)
  (void)MPI_Probe(source, tag, mpi_comm, &status.mpi_status);
    d220:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d224:	48 8d 48 08          	lea    0x8(%rax),%rcx
    d228:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    d22c:	48 8b 50 08          	mov    0x8(%rax),%rdx
    d230:	8b 75 f0             	mov    -0x10(%rbp),%esi
    d233:	8b 45 f4             	mov    -0xc(%rbp),%eax
    d236:	89 c7                	mov    %eax,%edi
    d238:	e8 43 cd ff ff       	call   9f80 <MPI_Probe@plt>
}
    d23d:	90                   	nop
    d23e:	c9                   	leave  
    d23f:	c3                   	ret    

000000000000d240 <_ZNK3MPI4Comm5ProbeEii>:

inline void
MPI::Comm::Probe(int source, int tag) const
{
    d240:	f3 0f 1e fa          	endbr64 
    d244:	55                   	push   %rbp
    d245:	48 89 e5             	mov    %rsp,%rbp
    d248:	48 83 ec 10          	sub    $0x10,%rsp
    d24c:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    d250:	89 75 f4             	mov    %esi,-0xc(%rbp)
    d253:	89 55 f0             	mov    %edx,-0x10(%rbp)
  (void)MPI_Probe(source, tag, mpi_comm, MPI_STATUS_IGNORE);
    d256:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    d25a:	48 8b 50 08          	mov    0x8(%rax),%rdx
    d25e:	8b 75 f0             	mov    -0x10(%rbp),%esi
    d261:	8b 45 f4             	mov    -0xc(%rbp),%eax
    d264:	b9 00 00 00 00       	mov    $0x0,%ecx
    d269:	89 c7                	mov    %eax,%edi
    d26b:	e8 10 cd ff ff       	call   9f80 <MPI_Probe@plt>
}
    d270:	90                   	nop
    d271:	c9                   	leave  
    d272:	c3                   	ret    
    d273:	90                   	nop

000000000000d274 <_ZNK3MPI4Comm9Send_initEPKviRKNS_8DatatypeEii>:

inline MPI::Prequest
MPI::Comm::Send_init(const void *buf, int count,
		     const MPI::Datatype & datatype, int dest, int tag) const
{
    d274:	f3 0f 1e fa          	endbr64 
    d278:	55                   	push   %rbp
    d279:	48 89 e5             	mov    %rsp,%rbp
    d27c:	53                   	push   %rbx
    d27d:	48 83 ec 48          	sub    $0x48,%rsp
    d281:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d285:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d289:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d28d:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d290:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d294:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d298:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d29f:	00 00 
    d2a1:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d2a5:	31 c0                	xor    %eax,%eax
  MPI_Request request;
  (void)MPI_Send_init(const_cast<void *>(buf), count, datatype,
		      dest, tag, mpi_comm, &request);
    d2a7:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d2ab:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Send_init(const_cast<void *>(buf), count, datatype,
    d2af:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d2b3:	48 89 c7             	mov    %rax,%rdi
    d2b6:	e8 c3 e5 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d2bb:	48 89 c7             	mov    %rax,%rdi
    d2be:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d2c1:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d2c4:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d2c8:	48 83 ec 08          	sub    $0x8,%rsp
    d2cc:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d2d0:	51                   	push   %rcx
    d2d1:	49 89 d9             	mov    %rbx,%r9
    d2d4:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d2d8:	89 d1                	mov    %edx,%ecx
    d2da:	48 89 fa             	mov    %rdi,%rdx
    d2dd:	48 89 c7             	mov    %rax,%rdi
    d2e0:	e8 bb d1 ff ff       	call   a4a0 <MPI_Send_init@plt>
    d2e5:	48 83 c4 10          	add    $0x10,%rsp
  return request;
    d2e9:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    d2ed:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d2f1:	48 89 d6             	mov    %rdx,%rsi
    d2f4:	48 89 c7             	mov    %rax,%rdi
    d2f7:	e8 98 e6 ff ff       	call   b994 <_ZN3MPI8PrequestC1ERKP14ompi_request_t>
}
    d2fc:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d300:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d307:	00 00 
    d309:	74 05                	je     d310 <_ZNK3MPI4Comm9Send_initEPKviRKNS_8DatatypeEii+0x9c>
    d30b:	e8 d0 d1 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    d310:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d314:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d318:	c9                   	leave  
    d319:	c3                   	ret    

000000000000d31a <_ZNK3MPI4Comm10Bsend_initEPKviRKNS_8DatatypeEii>:

inline MPI::Prequest
MPI::Comm::Bsend_init(const void *buf, int count,
		      const MPI::Datatype & datatype, int dest, int tag) const
{
    d31a:	f3 0f 1e fa          	endbr64 
    d31e:	55                   	push   %rbp
    d31f:	48 89 e5             	mov    %rsp,%rbp
    d322:	53                   	push   %rbx
    d323:	48 83 ec 48          	sub    $0x48,%rsp
    d327:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d32b:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d32f:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d333:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d336:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d33a:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d33e:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d345:	00 00 
    d347:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d34b:	31 c0                	xor    %eax,%eax
  MPI_Request request;
  (void)MPI_Bsend_init(const_cast<void *>(buf), count, datatype,
		       dest, tag, mpi_comm, &request);
    d34d:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d351:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Bsend_init(const_cast<void *>(buf), count, datatype,
    d355:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d359:	48 89 c7             	mov    %rax,%rdi
    d35c:	e8 1d e5 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d361:	48 89 c7             	mov    %rax,%rdi
    d364:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d367:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d36a:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d36e:	48 83 ec 08          	sub    $0x8,%rsp
    d372:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d376:	51                   	push   %rcx
    d377:	49 89 d9             	mov    %rbx,%r9
    d37a:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d37e:	89 d1                	mov    %edx,%ecx
    d380:	48 89 fa             	mov    %rdi,%rdx
    d383:	48 89 c7             	mov    %rax,%rdi
    d386:	e8 b5 ce ff ff       	call   a240 <MPI_Bsend_init@plt>
    d38b:	48 83 c4 10          	add    $0x10,%rsp
  return request;
    d38f:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    d393:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d397:	48 89 d6             	mov    %rdx,%rsi
    d39a:	48 89 c7             	mov    %rax,%rdi
    d39d:	e8 f2 e5 ff ff       	call   b994 <_ZN3MPI8PrequestC1ERKP14ompi_request_t>
}
    d3a2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d3a6:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d3ad:	00 00 
    d3af:	74 05                	je     d3b6 <_ZNK3MPI4Comm10Bsend_initEPKviRKNS_8DatatypeEii+0x9c>
    d3b1:	e8 2a d1 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    d3b6:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d3ba:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d3be:	c9                   	leave  
    d3bf:	c3                   	ret    

000000000000d3c0 <_ZNK3MPI4Comm10Ssend_initEPKviRKNS_8DatatypeEii>:

inline MPI::Prequest
MPI::Comm::Ssend_init(const void *buf, int count,
		      const MPI::Datatype & datatype, int dest, int tag) const
{
    d3c0:	f3 0f 1e fa          	endbr64 
    d3c4:	55                   	push   %rbp
    d3c5:	48 89 e5             	mov    %rsp,%rbp
    d3c8:	53                   	push   %rbx
    d3c9:	48 83 ec 48          	sub    $0x48,%rsp
    d3cd:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d3d1:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d3d5:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d3d9:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d3dc:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d3e0:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d3e4:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d3eb:	00 00 
    d3ed:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d3f1:	31 c0                	xor    %eax,%eax
  MPI_Request request;
  (void)MPI_Ssend_init(const_cast<void *>(buf), count, datatype,
		       dest, tag, mpi_comm, &request);
    d3f3:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d3f7:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Ssend_init(const_cast<void *>(buf), count, datatype,
    d3fb:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d3ff:	48 89 c7             	mov    %rax,%rdi
    d402:	e8 77 e4 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d407:	48 89 c7             	mov    %rax,%rdi
    d40a:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d40d:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d410:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d414:	48 83 ec 08          	sub    $0x8,%rsp
    d418:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d41c:	51                   	push   %rcx
    d41d:	49 89 d9             	mov    %rbx,%r9
    d420:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d424:	89 d1                	mov    %edx,%ecx
    d426:	48 89 fa             	mov    %rdi,%rdx
    d429:	48 89 c7             	mov    %rax,%rdi
    d42c:	e8 6f cd ff ff       	call   a1a0 <MPI_Ssend_init@plt>
    d431:	48 83 c4 10          	add    $0x10,%rsp
  return request;
    d435:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    d439:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d43d:	48 89 d6             	mov    %rdx,%rsi
    d440:	48 89 c7             	mov    %rax,%rdi
    d443:	e8 4c e5 ff ff       	call   b994 <_ZN3MPI8PrequestC1ERKP14ompi_request_t>
}
    d448:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d44c:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d453:	00 00 
    d455:	74 05                	je     d45c <_ZNK3MPI4Comm10Ssend_initEPKviRKNS_8DatatypeEii+0x9c>
    d457:	e8 84 d0 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    d45c:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d460:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d464:	c9                   	leave  
    d465:	c3                   	ret    

000000000000d466 <_ZNK3MPI4Comm10Rsend_initEPKviRKNS_8DatatypeEii>:

inline MPI::Prequest
MPI::Comm::Rsend_init(const void *buf, int count,
		      const MPI::Datatype & datatype, int dest, int tag) const
{
    d466:	f3 0f 1e fa          	endbr64 
    d46a:	55                   	push   %rbp
    d46b:	48 89 e5             	mov    %rsp,%rbp
    d46e:	53                   	push   %rbx
    d46f:	48 83 ec 48          	sub    $0x48,%rsp
    d473:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d477:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d47b:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d47f:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d482:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d486:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d48a:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d491:	00 00 
    d493:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d497:	31 c0                	xor    %eax,%eax
  MPI_Request request;
  (void)MPI_Rsend_init(const_cast<void *>(buf), count,  datatype,
		       dest, tag, mpi_comm, &request);
    d499:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d49d:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Rsend_init(const_cast<void *>(buf), count,  datatype,
    d4a1:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d4a5:	48 89 c7             	mov    %rax,%rdi
    d4a8:	e8 d1 e3 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d4ad:	48 89 c7             	mov    %rax,%rdi
    d4b0:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d4b3:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d4b6:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d4ba:	48 83 ec 08          	sub    $0x8,%rsp
    d4be:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d4c2:	51                   	push   %rcx
    d4c3:	49 89 d9             	mov    %rbx,%r9
    d4c6:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d4ca:	89 d1                	mov    %edx,%ecx
    d4cc:	48 89 fa             	mov    %rdi,%rdx
    d4cf:	48 89 c7             	mov    %rax,%rdi
    d4d2:	e8 49 cd ff ff       	call   a220 <MPI_Rsend_init@plt>
    d4d7:	48 83 c4 10          	add    $0x10,%rsp
  return request;
    d4db:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    d4df:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d4e3:	48 89 d6             	mov    %rdx,%rsi
    d4e6:	48 89 c7             	mov    %rax,%rdi
    d4e9:	e8 a6 e4 ff ff       	call   b994 <_ZN3MPI8PrequestC1ERKP14ompi_request_t>
}
    d4ee:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d4f2:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d4f9:	00 00 
    d4fb:	74 05                	je     d502 <_ZNK3MPI4Comm10Rsend_initEPKviRKNS_8DatatypeEii+0x9c>
    d4fd:	e8 de cf ff ff       	call   a4e0 <__stack_chk_fail@plt>
    d502:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d506:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d50a:	c9                   	leave  
    d50b:	c3                   	ret    

000000000000d50c <_ZNK3MPI4Comm9Recv_initEPviRKNS_8DatatypeEii>:

inline MPI::Prequest
MPI::Comm::Recv_init(void *buf, int count,
		     const MPI::Datatype & datatype, int source, int tag) const
{
    d50c:	f3 0f 1e fa          	endbr64 
    d510:	55                   	push   %rbp
    d511:	48 89 e5             	mov    %rsp,%rbp
    d514:	53                   	push   %rbx
    d515:	48 83 ec 48          	sub    $0x48,%rsp
    d519:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d51d:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d521:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d525:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d528:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d52c:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d530:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d537:	00 00 
    d539:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d53d:	31 c0                	xor    %eax,%eax
  MPI_Request request;
  (void)MPI_Recv_init(buf, count, datatype, source,
		      tag, mpi_comm, &request);
    d53f:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d543:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Recv_init(buf, count, datatype, source,
    d547:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d54b:	48 89 c7             	mov    %rax,%rdi
    d54e:	e8 2b e3 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d553:	48 89 c7             	mov    %rax,%rdi
    d556:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d559:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d55c:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d560:	48 83 ec 08          	sub    $0x8,%rsp
    d564:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d568:	51                   	push   %rcx
    d569:	49 89 d9             	mov    %rbx,%r9
    d56c:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d570:	89 d1                	mov    %edx,%ecx
    d572:	48 89 fa             	mov    %rdi,%rdx
    d575:	48 89 c7             	mov    %rax,%rdi
    d578:	e8 a3 c8 ff ff       	call   9e20 <MPI_Recv_init@plt>
    d57d:	48 83 c4 10          	add    $0x10,%rsp
  return request;
    d581:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    d585:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d589:	48 89 d6             	mov    %rdx,%rsi
    d58c:	48 89 c7             	mov    %rax,%rdi
    d58f:	e8 00 e4 ff ff       	call   b994 <_ZN3MPI8PrequestC1ERKP14ompi_request_t>
}
    d594:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d598:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d59f:	00 00 
    d5a1:	74 05                	je     d5a8 <_ZNK3MPI4Comm9Recv_initEPviRKNS_8DatatypeEii+0x9c>
    d5a3:	e8 38 cf ff ff       	call   a4e0 <__stack_chk_fail@plt>
    d5a8:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d5ac:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d5b0:	c9                   	leave  
    d5b1:	c3                   	ret    

000000000000d5b2 <_ZNK3MPI4Comm8SendrecvEPKviRKNS_8DatatypeEiiPviS5_iiRNS_6StatusE>:
MPI::Comm::Sendrecv(const void *sendbuf, int sendcount,
		    const MPI::Datatype & sendtype, int dest, int sendtag,
		    void *recvbuf, int recvcount,
		    const MPI::Datatype & recvtype, int source,
		    int recvtag, MPI::Status & status) const
{
    d5b2:	f3 0f 1e fa          	endbr64 
    d5b6:	55                   	push   %rbp
    d5b7:	48 89 e5             	mov    %rsp,%rbp
    d5ba:	41 55                	push   %r13
    d5bc:	41 54                	push   %r12
    d5be:	53                   	push   %rbx
    d5bf:	48 83 ec 38          	sub    $0x38,%rsp
    d5c3:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d5c7:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d5cb:	89 55 cc             	mov    %edx,-0x34(%rbp)
    d5ce:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    d5d2:	44 89 45 c8          	mov    %r8d,-0x38(%rbp)
    d5d6:	44 89 4d bc          	mov    %r9d,-0x44(%rbp)
  (void)MPI_Sendrecv(const_cast<void *>(sendbuf), sendcount,
    d5da:	48 8b 45 38          	mov    0x38(%rbp),%rax
    d5de:	4c 8d 68 08          	lea    0x8(%rax),%r13
		     sendtype,
		     dest, sendtag, recvbuf, recvcount,
		     recvtype,
		     source, recvtag, mpi_comm, &status.mpi_status);
    d5e2:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d5e6:	4c 8b 60 08          	mov    0x8(%rax),%r12
  (void)MPI_Sendrecv(const_cast<void *>(sendbuf), sendcount,
    d5ea:	48 8b 45 20          	mov    0x20(%rbp),%rax
    d5ee:	48 89 c7             	mov    %rax,%rdi
    d5f1:	e8 88 e2 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d5f6:	48 89 c3             	mov    %rax,%rbx
    d5f9:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    d5fd:	48 89 c7             	mov    %rax,%rdi
    d600:	e8 79 e2 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d605:	48 89 c7             	mov    %rax,%rdi
    d608:	44 8b 45 bc          	mov    -0x44(%rbp),%r8d
    d60c:	8b 55 c8             	mov    -0x38(%rbp),%edx
    d60f:	8b 75 cc             	mov    -0x34(%rbp),%esi
    d612:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d616:	41 55                	push   %r13
    d618:	41 54                	push   %r12
    d61a:	8b 4d 30             	mov    0x30(%rbp),%ecx
    d61d:	51                   	push   %rcx
    d61e:	8b 4d 28             	mov    0x28(%rbp),%ecx
    d621:	51                   	push   %rcx
    d622:	53                   	push   %rbx
    d623:	8b 4d 18             	mov    0x18(%rbp),%ecx
    d626:	51                   	push   %rcx
    d627:	4c 8b 4d 10          	mov    0x10(%rbp),%r9
    d62b:	89 d1                	mov    %edx,%ecx
    d62d:	48 89 fa             	mov    %rdi,%rdx
    d630:	48 89 c7             	mov    %rax,%rdi
    d633:	e8 28 cc ff ff       	call   a260 <MPI_Sendrecv@plt>
    d638:	48 83 c4 30          	add    $0x30,%rsp
}
    d63c:	90                   	nop
    d63d:	48 8d 65 e8          	lea    -0x18(%rbp),%rsp
    d641:	5b                   	pop    %rbx
    d642:	41 5c                	pop    %r12
    d644:	41 5d                	pop    %r13
    d646:	5d                   	pop    %rbp
    d647:	c3                   	ret    

000000000000d648 <_ZNK3MPI4Comm8SendrecvEPKviRKNS_8DatatypeEiiPviS5_ii>:
MPI::Comm::Sendrecv(const void *sendbuf, int sendcount,
		    const MPI::Datatype & sendtype, int dest, int sendtag,
		    void *recvbuf, int recvcount,
		    const MPI::Datatype & recvtype, int source,
		    int recvtag) const
{
    d648:	f3 0f 1e fa          	endbr64 
    d64c:	55                   	push   %rbp
    d64d:	48 89 e5             	mov    %rsp,%rbp
    d650:	41 54                	push   %r12
    d652:	53                   	push   %rbx
    d653:	48 83 ec 30          	sub    $0x30,%rsp
    d657:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d65b:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    d65f:	89 55 dc             	mov    %edx,-0x24(%rbp)
    d662:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    d666:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    d66a:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
  (void)MPI_Sendrecv(const_cast<void *>(sendbuf), sendcount,
		     sendtype,
		     dest, sendtag, recvbuf, recvcount,
		     recvtype,
		     source, recvtag, mpi_comm, MPI_STATUS_IGNORE);
    d66e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d672:	4c 8b 60 08          	mov    0x8(%rax),%r12
  (void)MPI_Sendrecv(const_cast<void *>(sendbuf), sendcount,
    d676:	48 8b 45 20          	mov    0x20(%rbp),%rax
    d67a:	48 89 c7             	mov    %rax,%rdi
    d67d:	e8 fc e1 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d682:	48 89 c3             	mov    %rax,%rbx
    d685:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d689:	48 89 c7             	mov    %rax,%rdi
    d68c:	e8 ed e1 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d691:	48 89 c7             	mov    %rax,%rdi
    d694:	44 8b 45 cc          	mov    -0x34(%rbp),%r8d
    d698:	8b 55 d8             	mov    -0x28(%rbp),%edx
    d69b:	8b 75 dc             	mov    -0x24(%rbp),%esi
    d69e:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    d6a2:	6a 00                	push   $0x0
    d6a4:	41 54                	push   %r12
    d6a6:	8b 4d 30             	mov    0x30(%rbp),%ecx
    d6a9:	51                   	push   %rcx
    d6aa:	8b 4d 28             	mov    0x28(%rbp),%ecx
    d6ad:	51                   	push   %rcx
    d6ae:	53                   	push   %rbx
    d6af:	8b 4d 18             	mov    0x18(%rbp),%ecx
    d6b2:	51                   	push   %rcx
    d6b3:	4c 8b 4d 10          	mov    0x10(%rbp),%r9
    d6b7:	89 d1                	mov    %edx,%ecx
    d6b9:	48 89 fa             	mov    %rdi,%rdx
    d6bc:	48 89 c7             	mov    %rax,%rdi
    d6bf:	e8 9c cb ff ff       	call   a260 <MPI_Sendrecv@plt>
    d6c4:	48 83 c4 30          	add    $0x30,%rsp
}
    d6c8:	90                   	nop
    d6c9:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    d6cd:	5b                   	pop    %rbx
    d6ce:	41 5c                	pop    %r12
    d6d0:	5d                   	pop    %rbp
    d6d1:	c3                   	ret    

000000000000d6d2 <_ZNK3MPI4Comm16Sendrecv_replaceEPviRKNS_8DatatypeEiiiiRNS_6StatusE>:
inline void
MPI::Comm::Sendrecv_replace(void *buf, int count,
			    const MPI::Datatype & datatype, int dest,
			    int sendtag, int source,
			    int recvtag, MPI::Status & status) const
{
    d6d2:	f3 0f 1e fa          	endbr64 
    d6d6:	55                   	push   %rbp
    d6d7:	48 89 e5             	mov    %rsp,%rbp
    d6da:	41 54                	push   %r12
    d6dc:	53                   	push   %rbx
    d6dd:	48 83 ec 30          	sub    $0x30,%rsp
    d6e1:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d6e5:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    d6e9:	89 55 dc             	mov    %edx,-0x24(%rbp)
    d6ec:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    d6f0:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    d6f4:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
  (void)MPI_Sendrecv_replace(buf, count, datatype, dest,
    d6f8:	48 8b 45 20          	mov    0x20(%rbp),%rax
    d6fc:	4c 8d 60 08          	lea    0x8(%rax),%r12
			     sendtag, source, recvtag, mpi_comm,
    d700:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d704:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Sendrecv_replace(buf, count, datatype, dest,
    d708:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d70c:	48 89 c7             	mov    %rax,%rdi
    d70f:	e8 6a e1 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d714:	48 89 c7             	mov    %rax,%rdi
    d717:	44 8b 45 cc          	mov    -0x34(%rbp),%r8d
    d71b:	8b 55 d8             	mov    -0x28(%rbp),%edx
    d71e:	8b 75 dc             	mov    -0x24(%rbp),%esi
    d721:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    d725:	48 83 ec 08          	sub    $0x8,%rsp
    d729:	41 54                	push   %r12
    d72b:	53                   	push   %rbx
    d72c:	8b 4d 18             	mov    0x18(%rbp),%ecx
    d72f:	51                   	push   %rcx
    d730:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
    d734:	89 d1                	mov    %edx,%ecx
    d736:	48 89 fa             	mov    %rdi,%rdx
    d739:	48 89 c7             	mov    %rax,%rdi
    d73c:	e8 1f cf ff ff       	call   a660 <MPI_Sendrecv_replace@plt>
    d741:	48 83 c4 20          	add    $0x20,%rsp
			     &status.mpi_status);
}
    d745:	90                   	nop
    d746:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    d74a:	5b                   	pop    %rbx
    d74b:	41 5c                	pop    %r12
    d74d:	5d                   	pop    %rbp
    d74e:	c3                   	ret    
    d74f:	90                   	nop

000000000000d750 <_ZNK3MPI4Comm16Sendrecv_replaceEPviRKNS_8DatatypeEiiii>:
inline void
MPI::Comm::Sendrecv_replace(void *buf, int count,
			    const MPI::Datatype & datatype, int dest,
			    int sendtag, int source,
			    int recvtag) const
{
    d750:	f3 0f 1e fa          	endbr64 
    d754:	55                   	push   %rbp
    d755:	48 89 e5             	mov    %rsp,%rbp
    d758:	53                   	push   %rbx
    d759:	48 83 ec 38          	sub    $0x38,%rsp
    d75d:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d761:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    d765:	89 55 dc             	mov    %edx,-0x24(%rbp)
    d768:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    d76c:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    d770:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
  (void)MPI_Sendrecv_replace(buf, count, datatype, dest,
			     sendtag, source, recvtag, mpi_comm,
    d774:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d778:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void)MPI_Sendrecv_replace(buf, count, datatype, dest,
    d77c:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d780:	48 89 c7             	mov    %rax,%rdi
    d783:	e8 f6 e0 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d788:	48 89 c7             	mov    %rax,%rdi
    d78b:	44 8b 45 cc          	mov    -0x34(%rbp),%r8d
    d78f:	8b 55 d8             	mov    -0x28(%rbp),%edx
    d792:	8b 75 dc             	mov    -0x24(%rbp),%esi
    d795:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    d799:	48 83 ec 08          	sub    $0x8,%rsp
    d79d:	6a 00                	push   $0x0
    d79f:	53                   	push   %rbx
    d7a0:	8b 4d 18             	mov    0x18(%rbp),%ecx
    d7a3:	51                   	push   %rcx
    d7a4:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
    d7a8:	89 d1                	mov    %edx,%ecx
    d7aa:	48 89 fa             	mov    %rdi,%rdx
    d7ad:	48 89 c7             	mov    %rax,%rdi
    d7b0:	e8 ab ce ff ff       	call   a660 <MPI_Sendrecv_replace@plt>
    d7b5:	48 83 c4 20          	add    $0x20,%rsp
			     MPI_STATUS_IGNORE);
}
    d7b9:	90                   	nop
    d7ba:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d7be:	c9                   	leave  
    d7bf:	c3                   	ret    

000000000000d7c0 <_ZNK3MPI4Comm9Get_groupEv>:
// Groups, Contexts, and Communicators
//

inline MPI::Group
MPI::Comm::Get_group() const
{
    d7c0:	f3 0f 1e fa          	endbr64 
    d7c4:	55                   	push   %rbp
    d7c5:	48 89 e5             	mov    %rsp,%rbp
    d7c8:	48 83 ec 20          	sub    $0x20,%rsp
    d7cc:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d7d0:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    d7d4:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d7db:	00 00 
    d7dd:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    d7e1:	31 c0                	xor    %eax,%eax
  MPI_Group group;
  (void)MPI_Comm_group(mpi_comm, &group);
    d7e3:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    d7e7:	48 8b 40 08          	mov    0x8(%rax),%rax
    d7eb:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    d7ef:	48 89 d6             	mov    %rdx,%rsi
    d7f2:	48 89 c7             	mov    %rax,%rdi
    d7f5:	e8 06 cb ff ff       	call   a300 <MPI_Comm_group@plt>
  return group;
    d7fa:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    d7fe:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d802:	48 89 d6             	mov    %rdx,%rsi
    d805:	48 89 c7             	mov    %rax,%rdi
    d808:	e8 7f e2 ff ff       	call   ba8c <_ZN3MPI5GroupC1EP12ompi_group_t>
}
    d80d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    d811:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d818:	00 00 
    d81a:	74 05                	je     d821 <_ZNK3MPI4Comm9Get_groupEv+0x61>
    d81c:	e8 bf cc ff ff       	call   a4e0 <__stack_chk_fail@plt>
    d821:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d825:	c9                   	leave  
    d826:	c3                   	ret    
    d827:	90                   	nop

000000000000d828 <_ZNK3MPI4Comm8Get_sizeEv>:

inline int
MPI::Comm::Get_size() const
{
    d828:	f3 0f 1e fa          	endbr64 
    d82c:	55                   	push   %rbp
    d82d:	48 89 e5             	mov    %rsp,%rbp
    d830:	48 83 ec 20          	sub    $0x20,%rsp
    d834:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d838:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d83f:	00 00 
    d841:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    d845:	31 c0                	xor    %eax,%eax
  int size;
  (void)MPI_Comm_size (mpi_comm, &size);
    d847:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d84b:	48 8b 40 08          	mov    0x8(%rax),%rax
    d84f:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    d853:	48 89 d6             	mov    %rdx,%rsi
    d856:	48 89 c7             	mov    %rax,%rdi
    d859:	e8 92 cc ff ff       	call   a4f0 <MPI_Comm_size@plt>
  return size;
    d85e:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    d861:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    d865:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    d86c:	00 00 
    d86e:	74 05                	je     d875 <_ZNK3MPI4Comm8Get_sizeEv+0x4d>
    d870:	e8 6b cc ff ff       	call   a4e0 <__stack_chk_fail@plt>
    d875:	c9                   	leave  
    d876:	c3                   	ret    
    d877:	90                   	nop

000000000000d878 <_ZNK3MPI4Comm8Get_rankEv>:

inline int
MPI::Comm::Get_rank() const
{
    d878:	f3 0f 1e fa          	endbr64 
    d87c:	55                   	push   %rbp
    d87d:	48 89 e5             	mov    %rsp,%rbp
    d880:	48 83 ec 20          	sub    $0x20,%rsp
    d884:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d888:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d88f:	00 00 
    d891:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    d895:	31 c0                	xor    %eax,%eax
  int myrank;
  (void)MPI_Comm_rank (mpi_comm, &myrank);
    d897:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d89b:	48 8b 40 08          	mov    0x8(%rax),%rax
    d89f:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    d8a3:	48 89 d6             	mov    %rdx,%rsi
    d8a6:	48 89 c7             	mov    %rax,%rdi
    d8a9:	e8 52 ce ff ff       	call   a700 <MPI_Comm_rank@plt>
  return myrank;
    d8ae:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    d8b1:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    d8b5:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    d8bc:	00 00 
    d8be:	74 05                	je     d8c5 <_ZNK3MPI4Comm8Get_rankEv+0x4d>
    d8c0:	e8 1b cc ff ff       	call   a4e0 <__stack_chk_fail@plt>
    d8c5:	c9                   	leave  
    d8c6:	c3                   	ret    
    d8c7:	90                   	nop

000000000000d8c8 <_ZN3MPI4Comm4FreeEv>:
  return result;
}

inline void
MPI::Comm::Free(void)
{
    d8c8:	f3 0f 1e fa          	endbr64 
    d8cc:	55                   	push   %rbp
    d8cd:	48 89 e5             	mov    %rsp,%rbp
    d8d0:	48 83 ec 10          	sub    $0x10,%rsp
    d8d4:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    (void)MPI_Comm_free(&mpi_comm);
    d8d8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    d8dc:	48 83 c0 08          	add    $0x8,%rax
    d8e0:	48 89 c7             	mov    %rax,%rdi
    d8e3:	e8 78 cf ff ff       	call   a860 <MPI_Comm_free@plt>
}
    d8e8:	90                   	nop
    d8e9:	c9                   	leave  
    d8ea:	c3                   	ret    
    d8eb:	90                   	nop

000000000000d8ec <_ZNK3MPI4Comm8Is_interEv>:

inline bool
MPI::Comm::Is_inter() const
{
    d8ec:	f3 0f 1e fa          	endbr64 
    d8f0:	55                   	push   %rbp
    d8f1:	48 89 e5             	mov    %rsp,%rbp
    d8f4:	48 83 ec 20          	sub    $0x20,%rsp
    d8f8:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d8fc:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d903:	00 00 
    d905:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    d909:	31 c0                	xor    %eax,%eax
  int t;
  (void)MPI_Comm_test_inter(mpi_comm, &t);
    d90b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d90f:	48 8b 40 08          	mov    0x8(%rax),%rax
    d913:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    d917:	48 89 d6             	mov    %rdx,%rsi
    d91a:	48 89 c7             	mov    %rax,%rdi
    d91d:	e8 fe d1 ff ff       	call   ab20 <MPI_Comm_test_inter@plt>
  return OPAL_INT_TO_BOOL(t);
    d922:	8b 45 f4             	mov    -0xc(%rbp),%eax
    d925:	85 c0                	test   %eax,%eax
    d927:	0f 95 c0             	setne  %al
}
    d92a:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    d92e:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    d935:	00 00 
    d937:	74 05                	je     d93e <_ZNK3MPI4Comm8Is_interEv+0x52>
    d939:	e8 a2 cb ff ff       	call   a4e0 <__stack_chk_fail@plt>
    d93e:	c9                   	leave  
    d93f:	c3                   	ret    

000000000000d940 <_ZNK3MPI4Comm7BarrierEv>:
// Collective Communication
//

inline void
MPI::Comm::Barrier() const
{
    d940:	f3 0f 1e fa          	endbr64 
    d944:	55                   	push   %rbp
    d945:	48 89 e5             	mov    %rsp,%rbp
    d948:	48 83 ec 10          	sub    $0x10,%rsp
    d94c:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    (void)MPI_Barrier(mpi_comm);
    d950:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    d954:	48 8b 40 08          	mov    0x8(%rax),%rax
    d958:	48 89 c7             	mov    %rax,%rdi
    d95b:	e8 80 c8 ff ff       	call   a1e0 <MPI_Barrier@plt>
}
    d960:	90                   	nop
    d961:	c9                   	leave  
    d962:	c3                   	ret    
    d963:	90                   	nop

000000000000d964 <_ZNK3MPI4Comm5BcastEPviRKNS_8DatatypeEi>:

inline void
MPI::Comm::Bcast(void *buffer, int count,
                 const MPI::Datatype& datatype, int root) const
{
    d964:	f3 0f 1e fa          	endbr64 
    d968:	55                   	push   %rbp
    d969:	48 89 e5             	mov    %rsp,%rbp
    d96c:	53                   	push   %rbx
    d96d:	48 83 ec 28          	sub    $0x28,%rsp
    d971:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d975:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    d979:	89 55 dc             	mov    %edx,-0x24(%rbp)
    d97c:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    d980:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    (void)MPI_Bcast(buffer, count, datatype, root, mpi_comm);
    d984:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d988:	48 8b 58 08          	mov    0x8(%rax),%rbx
    d98c:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d990:	48 89 c7             	mov    %rax,%rdi
    d993:	e8 e6 de ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d998:	48 89 c7             	mov    %rax,%rdi
    d99b:	8b 55 d8             	mov    -0x28(%rbp),%edx
    d99e:	8b 75 dc             	mov    -0x24(%rbp),%esi
    d9a1:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    d9a5:	49 89 d8             	mov    %rbx,%r8
    d9a8:	89 d1                	mov    %edx,%ecx
    d9aa:	48 89 fa             	mov    %rdi,%rdx
    d9ad:	48 89 c7             	mov    %rax,%rdi
    d9b0:	e8 bb c4 ff ff       	call   9e70 <MPI_Bcast@plt>
}
    d9b5:	90                   	nop
    d9b6:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d9ba:	c9                   	leave  
    d9bb:	c3                   	ret    

000000000000d9bc <_ZNK3MPI4Comm6GatherEPKviRKNS_8DatatypeEPviS5_i>:
inline void
MPI::Comm::Gather(const void *sendbuf, int sendcount,
                  const MPI::Datatype & sendtype,
                  void *recvbuf, int recvcount,
                  const MPI::Datatype & recvtype, int root) const
{
    d9bc:	f3 0f 1e fa          	endbr64 
    d9c0:	55                   	push   %rbp
    d9c1:	48 89 e5             	mov    %rsp,%rbp
    d9c4:	41 54                	push   %r12
    d9c6:	53                   	push   %rbx
    d9c7:	48 83 ec 30          	sub    $0x30,%rsp
    d9cb:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d9cf:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    d9d3:	89 55 dc             	mov    %edx,-0x24(%rbp)
    d9d6:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    d9da:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    d9de:	44 89 4d d8          	mov    %r9d,-0x28(%rbp)
    (void)MPI_Gather(const_cast<void *>(sendbuf), sendcount, sendtype,
                     recvbuf, recvcount, recvtype, root, mpi_comm);
    d9e2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d9e6:	48 8b 58 08          	mov    0x8(%rax),%rbx
    (void)MPI_Gather(const_cast<void *>(sendbuf), sendcount, sendtype,
    d9ea:	48 8b 7d 10          	mov    0x10(%rbp),%rdi
    d9ee:	e8 8b de ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d9f3:	49 89 c4             	mov    %rax,%r12
    d9f6:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d9fa:	48 89 c7             	mov    %rax,%rdi
    d9fd:	e8 7c de ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    da02:	48 89 c7             	mov    %rax,%rdi
    da05:	44 8b 45 d8          	mov    -0x28(%rbp),%r8d
    da09:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    da0d:	8b 75 dc             	mov    -0x24(%rbp),%esi
    da10:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    da14:	53                   	push   %rbx
    da15:	8b 4d 18             	mov    0x18(%rbp),%ecx
    da18:	51                   	push   %rcx
    da19:	4d 89 e1             	mov    %r12,%r9
    da1c:	48 89 d1             	mov    %rdx,%rcx
    da1f:	48 89 fa             	mov    %rdi,%rdx
    da22:	48 89 c7             	mov    %rax,%rdi
    da25:	e8 66 c8 ff ff       	call   a290 <MPI_Gather@plt>
    da2a:	48 83 c4 10          	add    $0x10,%rsp
}
    da2e:	90                   	nop
    da2f:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    da33:	5b                   	pop    %rbx
    da34:	41 5c                	pop    %r12
    da36:	5d                   	pop    %rbp
    da37:	c3                   	ret    

000000000000da38 <_ZNK3MPI4Comm7GathervEPKviRKNS_8DatatypeEPvPKiS8_S5_i>:
inline void
MPI::Comm::Gatherv(const void *sendbuf, int sendcount,
                   const MPI::Datatype & sendtype, void *recvbuf,
                   const int recvcounts[], const int displs[],
                   const MPI::Datatype & recvtype, int root) const
{
    da38:	f3 0f 1e fa          	endbr64 
    da3c:	55                   	push   %rbp
    da3d:	48 89 e5             	mov    %rsp,%rbp
    da40:	41 54                	push   %r12
    da42:	53                   	push   %rbx
    da43:	48 83 ec 30          	sub    $0x30,%rsp
    da47:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    da4b:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    da4f:	89 55 dc             	mov    %edx,-0x24(%rbp)
    da52:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    da56:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    da5a:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    (void)MPI_Gatherv(const_cast<void *>(sendbuf), sendcount,  sendtype,
                      recvbuf, const_cast<int *>(recvcounts),
                      const_cast<int *>(displs),
                      recvtype, root, mpi_comm);
    da5e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    da62:	4c 8b 60 08          	mov    0x8(%rax),%r12
    (void)MPI_Gatherv(const_cast<void *>(sendbuf), sendcount,  sendtype,
    da66:	48 8b 45 18          	mov    0x18(%rbp),%rax
    da6a:	48 89 c7             	mov    %rax,%rdi
    da6d:	e8 0c de ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    da72:	48 89 c3             	mov    %rax,%rbx
    da75:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    da79:	48 89 c7             	mov    %rax,%rdi
    da7c:	e8 fd dd ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    da81:	48 89 c7             	mov    %rax,%rdi
    da84:	4c 8b 45 c0          	mov    -0x40(%rbp),%r8
    da88:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    da8c:	8b 75 dc             	mov    -0x24(%rbp),%esi
    da8f:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    da93:	48 83 ec 08          	sub    $0x8,%rsp
    da97:	41 54                	push   %r12
    da99:	8b 4d 20             	mov    0x20(%rbp),%ecx
    da9c:	51                   	push   %rcx
    da9d:	53                   	push   %rbx
    da9e:	4c 8b 4d 10          	mov    0x10(%rbp),%r9
    daa2:	48 89 d1             	mov    %rdx,%rcx
    daa5:	48 89 fa             	mov    %rdi,%rdx
    daa8:	48 89 c7             	mov    %rax,%rdi
    daab:	e8 70 cf ff ff       	call   aa20 <MPI_Gatherv@plt>
    dab0:	48 83 c4 20          	add    $0x20,%rsp
}
    dab4:	90                   	nop
    dab5:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    dab9:	5b                   	pop    %rbx
    daba:	41 5c                	pop    %r12
    dabc:	5d                   	pop    %rbp
    dabd:	c3                   	ret    

000000000000dabe <_ZNK3MPI4Comm7ScatterEPKviRKNS_8DatatypeEPviS5_i>:
inline void
MPI::Comm::Scatter(const void *sendbuf, int sendcount,
                   const MPI::Datatype & sendtype,
                   void *recvbuf, int recvcount,
                   const MPI::Datatype & recvtype, int root) const
{
    dabe:	f3 0f 1e fa          	endbr64 
    dac2:	55                   	push   %rbp
    dac3:	48 89 e5             	mov    %rsp,%rbp
    dac6:	41 54                	push   %r12
    dac8:	53                   	push   %rbx
    dac9:	48 83 ec 30          	sub    $0x30,%rsp
    dacd:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    dad1:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    dad5:	89 55 dc             	mov    %edx,-0x24(%rbp)
    dad8:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    dadc:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    dae0:	44 89 4d d8          	mov    %r9d,-0x28(%rbp)
    (void)MPI_Scatter(const_cast<void *>(sendbuf), sendcount, sendtype,
                      recvbuf, recvcount, recvtype, root, mpi_comm);
    dae4:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    dae8:	48 8b 58 08          	mov    0x8(%rax),%rbx
    (void)MPI_Scatter(const_cast<void *>(sendbuf), sendcount, sendtype,
    daec:	48 8b 7d 10          	mov    0x10(%rbp),%rdi
    daf0:	e8 89 dd ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    daf5:	49 89 c4             	mov    %rax,%r12
    daf8:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    dafc:	48 89 c7             	mov    %rax,%rdi
    daff:	e8 7a dd ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    db04:	48 89 c7             	mov    %rax,%rdi
    db07:	44 8b 45 d8          	mov    -0x28(%rbp),%r8d
    db0b:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    db0f:	8b 75 dc             	mov    -0x24(%rbp),%esi
    db12:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    db16:	53                   	push   %rbx
    db17:	8b 4d 18             	mov    0x18(%rbp),%ecx
    db1a:	51                   	push   %rcx
    db1b:	4d 89 e1             	mov    %r12,%r9
    db1e:	48 89 d1             	mov    %rdx,%rcx
    db21:	48 89 fa             	mov    %rdi,%rdx
    db24:	48 89 c7             	mov    %rax,%rdi
    db27:	e8 44 c6 ff ff       	call   a170 <MPI_Scatter@plt>
    db2c:	48 83 c4 10          	add    $0x10,%rsp
}
    db30:	90                   	nop
    db31:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    db35:	5b                   	pop    %rbx
    db36:	41 5c                	pop    %r12
    db38:	5d                   	pop    %rbp
    db39:	c3                   	ret    

000000000000db3a <_ZNK3MPI4Comm8ScattervEPKvPKiS4_RKNS_8DatatypeEPviS7_i>:
inline void
MPI::Comm::Scatterv(const void *sendbuf, const int sendcounts[],
                    const int displs[], const MPI::Datatype & sendtype,
                    void *recvbuf, int recvcount,
                    const MPI::Datatype & recvtype, int root) const
{
    db3a:	f3 0f 1e fa          	endbr64 
    db3e:	55                   	push   %rbp
    db3f:	48 89 e5             	mov    %rsp,%rbp
    db42:	41 54                	push   %r12
    db44:	53                   	push   %rbx
    db45:	48 83 ec 30          	sub    $0x30,%rsp
    db49:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    db4d:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    db51:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    db55:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    db59:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    db5d:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    (void)MPI_Scatterv(const_cast<void *>(sendbuf),
                       const_cast<int *>(sendcounts),
                       const_cast<int *>(displs), sendtype,
                       recvbuf, recvcount, recvtype,
                       root, mpi_comm);
    db61:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    db65:	4c 8b 60 08          	mov    0x8(%rax),%r12
    (void)MPI_Scatterv(const_cast<void *>(sendbuf),
    db69:	48 8b 45 18          	mov    0x18(%rbp),%rax
    db6d:	48 89 c7             	mov    %rax,%rdi
    db70:	e8 09 dd ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    db75:	48 89 c3             	mov    %rax,%rbx
    db78:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    db7c:	48 89 c7             	mov    %rax,%rdi
    db7f:	e8 fa dc ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    db84:	48 89 c7             	mov    %rax,%rdi
    db87:	4c 8b 45 c0          	mov    -0x40(%rbp),%r8
    db8b:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    db8f:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    db93:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    db97:	48 83 ec 08          	sub    $0x8,%rsp
    db9b:	41 54                	push   %r12
    db9d:	8b 4d 20             	mov    0x20(%rbp),%ecx
    dba0:	51                   	push   %rcx
    dba1:	53                   	push   %rbx
    dba2:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
    dba6:	48 89 f9             	mov    %rdi,%rcx
    dba9:	48 89 c7             	mov    %rax,%rdi
    dbac:	e8 df c4 ff ff       	call   a090 <MPI_Scatterv@plt>
    dbb1:	48 83 c4 20          	add    $0x20,%rsp
}
    dbb5:	90                   	nop
    dbb6:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    dbba:	5b                   	pop    %rbx
    dbbb:	41 5c                	pop    %r12
    dbbd:	5d                   	pop    %rbp
    dbbe:	c3                   	ret    
    dbbf:	90                   	nop

000000000000dbc0 <_ZNK3MPI4Comm9AllgatherEPKviRKNS_8DatatypeEPviS5_>:

inline void
MPI::Comm::Allgather(const void *sendbuf, int sendcount,
                     const MPI::Datatype & sendtype, void *recvbuf,
                     int recvcount, const MPI::Datatype & recvtype) const
{
    dbc0:	f3 0f 1e fa          	endbr64 
    dbc4:	55                   	push   %rbp
    dbc5:	48 89 e5             	mov    %rsp,%rbp
    dbc8:	41 54                	push   %r12
    dbca:	53                   	push   %rbx
    dbcb:	48 83 ec 30          	sub    $0x30,%rsp
    dbcf:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    dbd3:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    dbd7:	89 55 dc             	mov    %edx,-0x24(%rbp)
    dbda:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    dbde:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    dbe2:	44 89 4d d8          	mov    %r9d,-0x28(%rbp)
    (void)MPI_Allgather(const_cast<void *>(sendbuf), sendcount,
                        sendtype, recvbuf, recvcount,
                        recvtype, mpi_comm);
    dbe6:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    dbea:	48 8b 58 08          	mov    0x8(%rax),%rbx
    (void)MPI_Allgather(const_cast<void *>(sendbuf), sendcount,
    dbee:	48 8b 7d 10          	mov    0x10(%rbp),%rdi
    dbf2:	e8 87 dc ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dbf7:	49 89 c4             	mov    %rax,%r12
    dbfa:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    dbfe:	48 89 c7             	mov    %rax,%rdi
    dc01:	e8 78 dc ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dc06:	48 89 c7             	mov    %rax,%rdi
    dc09:	8b 4d d8             	mov    -0x28(%rbp),%ecx
    dc0c:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    dc10:	8b 75 dc             	mov    -0x24(%rbp),%esi
    dc13:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    dc17:	48 83 ec 08          	sub    $0x8,%rsp
    dc1b:	53                   	push   %rbx
    dc1c:	4d 89 e1             	mov    %r12,%r9
    dc1f:	41 89 c8             	mov    %ecx,%r8d
    dc22:	48 89 d1             	mov    %rdx,%rcx
    dc25:	48 89 fa             	mov    %rdi,%rdx
    dc28:	48 89 c7             	mov    %rax,%rdi
    dc2b:	e8 e0 c1 ff ff       	call   9e10 <MPI_Allgather@plt>
    dc30:	48 83 c4 10          	add    $0x10,%rsp
}
    dc34:	90                   	nop
    dc35:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    dc39:	5b                   	pop    %rbx
    dc3a:	41 5c                	pop    %r12
    dc3c:	5d                   	pop    %rbp
    dc3d:	c3                   	ret    

000000000000dc3e <_ZNK3MPI4Comm10AllgathervEPKviRKNS_8DatatypeEPvPKiS8_S5_>:
inline void
MPI::Comm::Allgatherv(const void *sendbuf, int sendcount,
                      const MPI::Datatype & sendtype, void *recvbuf,
                      const int recvcounts[], const int displs[],
                      const MPI::Datatype & recvtype) const
{
    dc3e:	f3 0f 1e fa          	endbr64 
    dc42:	55                   	push   %rbp
    dc43:	48 89 e5             	mov    %rsp,%rbp
    dc46:	41 54                	push   %r12
    dc48:	53                   	push   %rbx
    dc49:	48 83 ec 30          	sub    $0x30,%rsp
    dc4d:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    dc51:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    dc55:	89 55 dc             	mov    %edx,-0x24(%rbp)
    dc58:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    dc5c:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    dc60:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    (void)MPI_Allgatherv(const_cast<void *>(sendbuf), sendcount,
                         sendtype, recvbuf,
                         const_cast<int *>(recvcounts),
                         const_cast<int *>(displs),
                         recvtype, mpi_comm);
    dc64:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    dc68:	4c 8b 60 08          	mov    0x8(%rax),%r12
    (void)MPI_Allgatherv(const_cast<void *>(sendbuf), sendcount,
    dc6c:	48 8b 45 18          	mov    0x18(%rbp),%rax
    dc70:	48 89 c7             	mov    %rax,%rdi
    dc73:	e8 06 dc ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dc78:	48 89 c3             	mov    %rax,%rbx
    dc7b:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    dc7f:	48 89 c7             	mov    %rax,%rdi
    dc82:	e8 f7 db ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dc87:	48 89 c7             	mov    %rax,%rdi
    dc8a:	48 8b 4d c0          	mov    -0x40(%rbp),%rcx
    dc8e:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    dc92:	8b 75 dc             	mov    -0x24(%rbp),%esi
    dc95:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    dc99:	41 54                	push   %r12
    dc9b:	53                   	push   %rbx
    dc9c:	4c 8b 4d 10          	mov    0x10(%rbp),%r9
    dca0:	49 89 c8             	mov    %rcx,%r8
    dca3:	48 89 d1             	mov    %rdx,%rcx
    dca6:	48 89 fa             	mov    %rdi,%rdx
    dca9:	48 89 c7             	mov    %rax,%rdi
    dcac:	e8 bf c8 ff ff       	call   a570 <MPI_Allgatherv@plt>
    dcb1:	48 83 c4 10          	add    $0x10,%rsp
}
    dcb5:	90                   	nop
    dcb6:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    dcba:	5b                   	pop    %rbx
    dcbb:	41 5c                	pop    %r12
    dcbd:	5d                   	pop    %rbp
    dcbe:	c3                   	ret    
    dcbf:	90                   	nop

000000000000dcc0 <_ZNK3MPI4Comm8AlltoallEPKviRKNS_8DatatypeEPviS5_>:

inline void
MPI::Comm::Alltoall(const void *sendbuf, int sendcount,
                    const MPI::Datatype & sendtype, void *recvbuf,
                    int recvcount, const MPI::Datatype & recvtype) const
{
    dcc0:	f3 0f 1e fa          	endbr64 
    dcc4:	55                   	push   %rbp
    dcc5:	48 89 e5             	mov    %rsp,%rbp
    dcc8:	41 54                	push   %r12
    dcca:	53                   	push   %rbx
    dccb:	48 83 ec 30          	sub    $0x30,%rsp
    dccf:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    dcd3:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    dcd7:	89 55 dc             	mov    %edx,-0x24(%rbp)
    dcda:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    dcde:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    dce2:	44 89 4d d8          	mov    %r9d,-0x28(%rbp)
    (void)MPI_Alltoall(const_cast<void *>(sendbuf), sendcount,
                       sendtype, recvbuf, recvcount,
                       recvtype, mpi_comm);
    dce6:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    dcea:	48 8b 58 08          	mov    0x8(%rax),%rbx
    (void)MPI_Alltoall(const_cast<void *>(sendbuf), sendcount,
    dcee:	48 8b 7d 10          	mov    0x10(%rbp),%rdi
    dcf2:	e8 87 db ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dcf7:	49 89 c4             	mov    %rax,%r12
    dcfa:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    dcfe:	48 89 c7             	mov    %rax,%rdi
    dd01:	e8 78 db ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dd06:	48 89 c7             	mov    %rax,%rdi
    dd09:	8b 4d d8             	mov    -0x28(%rbp),%ecx
    dd0c:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    dd10:	8b 75 dc             	mov    -0x24(%rbp),%esi
    dd13:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    dd17:	48 83 ec 08          	sub    $0x8,%rsp
    dd1b:	53                   	push   %rbx
    dd1c:	4d 89 e1             	mov    %r12,%r9
    dd1f:	41 89 c8             	mov    %ecx,%r8d
    dd22:	48 89 d1             	mov    %rdx,%rcx
    dd25:	48 89 fa             	mov    %rdi,%rdx
    dd28:	48 89 c7             	mov    %rax,%rdi
    dd2b:	e8 30 c3 ff ff       	call   a060 <MPI_Alltoall@plt>
    dd30:	48 83 c4 10          	add    $0x10,%rsp
}
    dd34:	90                   	nop
    dd35:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    dd39:	5b                   	pop    %rbx
    dd3a:	41 5c                	pop    %r12
    dd3c:	5d                   	pop    %rbp
    dd3d:	c3                   	ret    

000000000000dd3e <_ZNK3MPI4Comm9AlltoallvEPKvPKiS4_RKNS_8DatatypeEPvS4_S4_S7_>:
MPI::Comm::Alltoallv(const void *sendbuf, const int sendcounts[],
                     const int sdispls[], const MPI::Datatype & sendtype,
                     void *recvbuf, const int recvcounts[],
                     const int rdispls[],
                     const MPI::Datatype & recvtype) const
{
    dd3e:	f3 0f 1e fa          	endbr64 
    dd42:	55                   	push   %rbp
    dd43:	48 89 e5             	mov    %rsp,%rbp
    dd46:	41 54                	push   %r12
    dd48:	53                   	push   %rbx
    dd49:	48 83 ec 30          	sub    $0x30,%rsp
    dd4d:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    dd51:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    dd55:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    dd59:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    dd5d:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    dd61:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    (void)MPI_Alltoallv(const_cast<void *>(sendbuf),
                        const_cast<int *>(sendcounts),
			const_cast<int *>(sdispls), sendtype, recvbuf,
			const_cast<int *>(recvcounts),
                        const_cast<int *>(rdispls),
			recvtype,mpi_comm);
    dd65:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    dd69:	4c 8b 60 08          	mov    0x8(%rax),%r12
    (void)MPI_Alltoallv(const_cast<void *>(sendbuf),
    dd6d:	48 8b 45 20          	mov    0x20(%rbp),%rax
    dd71:	48 89 c7             	mov    %rax,%rdi
    dd74:	e8 05 db ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dd79:	48 89 c3             	mov    %rax,%rbx
    dd7c:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    dd80:	48 89 c7             	mov    %rax,%rdi
    dd83:	e8 f6 da ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dd88:	48 89 c1             	mov    %rax,%rcx
    dd8b:	48 8b 7d c0          	mov    -0x40(%rbp),%rdi
    dd8f:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    dd93:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    dd97:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    dd9b:	48 83 ec 08          	sub    $0x8,%rsp
    dd9f:	41 54                	push   %r12
    dda1:	53                   	push   %rbx
    dda2:	ff 75 18             	push   0x18(%rbp)
    dda5:	4c 8b 4d 10          	mov    0x10(%rbp),%r9
    dda9:	49 89 f8             	mov    %rdi,%r8
    ddac:	48 89 c7             	mov    %rax,%rdi
    ddaf:	e8 5c ca ff ff       	call   a810 <MPI_Alltoallv@plt>
    ddb4:	48 83 c4 20          	add    $0x20,%rsp
}
    ddb8:	90                   	nop
    ddb9:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    ddbd:	5b                   	pop    %rbx
    ddbe:	41 5c                	pop    %r12
    ddc0:	5d                   	pop    %rbp
    ddc1:	c3                   	ret    

000000000000ddc2 <_ZNK3MPI4Comm9AlltoallwEPKvPKiS4_PKNS_8DatatypeEPvS4_S4_S7_>:
MPI::Comm::Alltoallw(const void *sendbuf, const int sendcounts[],
                     const int sdispls[], const MPI::Datatype sendtypes[],
                     void *recvbuf, const int recvcounts[],
                     const int rdispls[],
                     const MPI::Datatype recvtypes[]) const
{
    ddc2:	f3 0f 1e fa          	endbr64 
    ddc6:	55                   	push   %rbp
    ddc7:	48 89 e5             	mov    %rsp,%rbp
    ddca:	53                   	push   %rbx
    ddcb:	48 83 ec 48          	sub    $0x48,%rsp
    ddcf:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    ddd3:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    ddd7:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    dddb:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    dddf:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    dde3:	4c 89 4d b0          	mov    %r9,-0x50(%rbp)
    const int comm_size = Get_size();
    dde7:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ddeb:	48 8b 00             	mov    (%rax),%rax
    ddee:	48 05 d8 00 00 00    	add    $0xd8,%rax
    ddf4:	48 8b 10             	mov    (%rax),%rdx
    ddf7:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ddfb:	48 89 c7             	mov    %rax,%rdi
    ddfe:	ff d2                	call   *%rdx
    de00:	89 45 e4             	mov    %eax,-0x1c(%rbp)
    MPI_Datatype *const data_type_tbl = new MPI_Datatype [2*comm_size];
    de03:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    de06:	01 c0                	add    %eax,%eax
    de08:	48 98                	cltq   
    de0a:	48 ba ff ff ff ff ff 	movabs $0xfffffffffffffff,%rdx
    de11:	ff ff 0f 
    de14:	48 39 d0             	cmp    %rdx,%rax
    de17:	77 20                	ja     de39 <_ZNK3MPI4Comm9AlltoallwEPKvPKiS4_PKNS_8DatatypeEPvS4_S4_S7_+0x77>
    de19:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    de1c:	01 c0                	add    %eax,%eax
    de1e:	48 98                	cltq   
    de20:	48 c1 e0 03          	shl    $0x3,%rax
    de24:	48 89 c7             	mov    %rax,%rdi
    de27:	e8 94 bf ff ff       	call   9dc0 <_Znam@plt>
    de2c:	48 89 45 e8          	mov    %rax,-0x18(%rbp)

    // This must be done because MPI::Datatype arrays cannot be
    // converted directly into MPI_Datatype arrays.
    for (int i_rank=0; i_rank < comm_size; i_rank++) {
    de30:	c7 45 e0 00 00 00 00 	movl   $0x0,-0x20(%rbp)
    de37:	eb 76                	jmp    deaf <_ZNK3MPI4Comm9AlltoallwEPKvPKiS4_PKNS_8DatatypeEPvS4_S4_S7_+0xed>
    MPI_Datatype *const data_type_tbl = new MPI_Datatype [2*comm_size];
    de39:	e8 b2 c8 ff ff       	call   a6f0 <__cxa_throw_bad_array_new_length@plt>
        data_type_tbl[i_rank] = sendtypes[i_rank];
    de3e:	8b 45 e0             	mov    -0x20(%rbp),%eax
    de41:	48 98                	cltq   
    de43:	48 c1 e0 04          	shl    $0x4,%rax
    de47:	48 89 c2             	mov    %rax,%rdx
    de4a:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    de4e:	48 01 d0             	add    %rdx,%rax
    de51:	8b 55 e0             	mov    -0x20(%rbp),%edx
    de54:	48 63 d2             	movslq %edx,%rdx
    de57:	48 8d 0c d5 00 00 00 	lea    0x0(,%rdx,8),%rcx
    de5e:	00 
    de5f:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    de63:	48 8d 1c 11          	lea    (%rcx,%rdx,1),%rbx
    de67:	48 89 c7             	mov    %rax,%rdi
    de6a:	e8 0f da ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    de6f:	48 89 03             	mov    %rax,(%rbx)
        data_type_tbl[i_rank + comm_size] = recvtypes[i_rank];
    de72:	8b 45 e0             	mov    -0x20(%rbp),%eax
    de75:	48 98                	cltq   
    de77:	48 c1 e0 04          	shl    $0x4,%rax
    de7b:	48 89 c2             	mov    %rax,%rdx
    de7e:	48 8b 45 20          	mov    0x20(%rbp),%rax
    de82:	48 01 d0             	add    %rdx,%rax
    de85:	8b 4d e0             	mov    -0x20(%rbp),%ecx
    de88:	8b 55 e4             	mov    -0x1c(%rbp),%edx
    de8b:	01 ca                	add    %ecx,%edx
    de8d:	48 63 d2             	movslq %edx,%rdx
    de90:	48 8d 0c d5 00 00 00 	lea    0x0(,%rdx,8),%rcx
    de97:	00 
    de98:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    de9c:	48 8d 1c 11          	lea    (%rcx,%rdx,1),%rbx
    dea0:	48 89 c7             	mov    %rax,%rdi
    dea3:	e8 d6 d9 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dea8:	48 89 03             	mov    %rax,(%rbx)
    for (int i_rank=0; i_rank < comm_size; i_rank++) {
    deab:	83 45 e0 01          	addl   $0x1,-0x20(%rbp)
    deaf:	8b 45 e0             	mov    -0x20(%rbp),%eax
    deb2:	3b 45 e4             	cmp    -0x1c(%rbp),%eax
    deb5:	7c 87                	jl     de3e <_ZNK3MPI4Comm9AlltoallwEPKvPKiS4_PKNS_8DatatypeEPvS4_S4_S7_+0x7c>
                        const_cast<int *>(sendcounts),
                        const_cast<int *>(sdispls),
                        data_type_tbl, recvbuf,
                        const_cast<int *>(recvcounts),
                        const_cast<int *>(rdispls),
                        &data_type_tbl[comm_size], mpi_comm);
    deb7:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    debb:	48 8b 78 08          	mov    0x8(%rax),%rdi
    debf:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    dec2:	48 98                	cltq   
    dec4:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
    decb:	00 
    (void)MPI_Alltoallw(const_cast<void *>(sendbuf),
    decc:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ded0:	4c 8d 04 02          	lea    (%rdx,%rax,1),%r8
    ded4:	4c 8b 55 b0          	mov    -0x50(%rbp),%r10
    ded8:	48 8b 4d e8          	mov    -0x18(%rbp),%rcx
    dedc:	48 8b 55 c0          	mov    -0x40(%rbp),%rdx
    dee0:	48 8b 75 c8          	mov    -0x38(%rbp),%rsi
    dee4:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    dee8:	48 83 ec 08          	sub    $0x8,%rsp
    deec:	57                   	push   %rdi
    deed:	41 50                	push   %r8
    deef:	ff 75 18             	push   0x18(%rbp)
    def2:	4c 8b 4d 10          	mov    0x10(%rbp),%r9
    def6:	4d 89 d0             	mov    %r10,%r8
    def9:	48 89 c7             	mov    %rax,%rdi
    defc:	e8 8f c7 ff ff       	call   a690 <MPI_Alltoallw@plt>
    df01:	48 83 c4 20          	add    $0x20,%rsp

    delete[] data_type_tbl;
    df05:	48 83 7d e8 00       	cmpq   $0x0,-0x18(%rbp)
    df0a:	74 0c                	je     df18 <_ZNK3MPI4Comm9AlltoallwEPKvPKiS4_PKNS_8DatatypeEPvS4_S4_S7_+0x156>
    df0c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    df10:	48 89 c7             	mov    %rax,%rdi
    df13:	e8 b8 c6 ff ff       	call   a5d0 <_ZdaPv@plt>
}
    df18:	90                   	nop
    df19:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    df1d:	c9                   	leave  
    df1e:	c3                   	ret    
    df1f:	90                   	nop

000000000000df20 <_ZNK3MPI4Comm6ReduceEPKvPviRKNS_8DatatypeERKNS_2OpEi>:

inline void
MPI::Comm::Reduce(const void *sendbuf, void *recvbuf, int count,
                  const MPI::Datatype & datatype, const MPI::Op& op,
                  int root) const
{
    df20:	f3 0f 1e fa          	endbr64 
    df24:	55                   	push   %rbp
    df25:	48 89 e5             	mov    %rsp,%rbp
    df28:	41 54                	push   %r12
    df2a:	53                   	push   %rbx
    df2b:	48 83 ec 30          	sub    $0x30,%rsp
    df2f:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    df33:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    df37:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    df3b:	89 4d d4             	mov    %ecx,-0x2c(%rbp)
    df3e:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    df42:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    (void)MPI_Reduce(const_cast<void *>(sendbuf), recvbuf, count, datatype, op, root, mpi_comm);
    df46:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    df4a:	48 8b 58 08          	mov    0x8(%rax),%rbx
    df4e:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    df52:	48 89 c7             	mov    %rax,%rdi
    df55:	e8 9a 1b 00 00       	call   faf4 <_ZNK3MPI2OpcvP9ompi_op_tEv>
    df5a:	49 89 c4             	mov    %rax,%r12
    df5d:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    df61:	48 89 c7             	mov    %rax,%rdi
    df64:	e8 15 d9 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    df69:	48 89 c1             	mov    %rax,%rcx
    df6c:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    df6f:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    df73:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    df77:	48 83 ec 08          	sub    $0x8,%rsp
    df7b:	53                   	push   %rbx
    df7c:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
    df80:	4d 89 e0             	mov    %r12,%r8
    df83:	48 89 c7             	mov    %rax,%rdi
    df86:	e8 85 c6 ff ff       	call   a610 <MPI_Reduce@plt>
    df8b:	48 83 c4 10          	add    $0x10,%rsp
}
    df8f:	90                   	nop
    df90:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    df94:	5b                   	pop    %rbx
    df95:	41 5c                	pop    %r12
    df97:	5d                   	pop    %rbp
    df98:	c3                   	ret    
    df99:	90                   	nop

000000000000df9a <_ZNK3MPI4Comm9AllreduceEPKvPviRKNS_8DatatypeERKNS_2OpE>:

inline void
MPI::Comm::Allreduce(const void *sendbuf, void *recvbuf, int count,
                     const MPI::Datatype & datatype, const MPI::Op& op) const
{
    df9a:	f3 0f 1e fa          	endbr64 
    df9e:	55                   	push   %rbp
    df9f:	48 89 e5             	mov    %rsp,%rbp
    dfa2:	41 54                	push   %r12
    dfa4:	53                   	push   %rbx
    dfa5:	48 83 ec 30          	sub    $0x30,%rsp
    dfa9:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    dfad:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    dfb1:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    dfb5:	89 4d d4             	mov    %ecx,-0x2c(%rbp)
    dfb8:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    dfbc:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    (void)MPI_Allreduce (const_cast<void *>(sendbuf), recvbuf, count, datatype,  op, mpi_comm);
    dfc0:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    dfc4:	4c 8b 60 08          	mov    0x8(%rax),%r12
    dfc8:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    dfcc:	48 89 c7             	mov    %rax,%rdi
    dfcf:	e8 20 1b 00 00       	call   faf4 <_ZNK3MPI2OpcvP9ompi_op_tEv>
    dfd4:	48 89 c3             	mov    %rax,%rbx
    dfd7:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    dfdb:	48 89 c7             	mov    %rax,%rdi
    dfde:	e8 9b d8 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dfe3:	48 89 c1             	mov    %rax,%rcx
    dfe6:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    dfe9:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    dfed:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    dff1:	4d 89 e1             	mov    %r12,%r9
    dff4:	49 89 d8             	mov    %rbx,%r8
    dff7:	48 89 c7             	mov    %rax,%rdi
    dffa:	e8 61 c5 ff ff       	call   a560 <MPI_Allreduce@plt>
}
    dfff:	90                   	nop
    e000:	48 83 c4 30          	add    $0x30,%rsp
    e004:	5b                   	pop    %rbx
    e005:	41 5c                	pop    %r12
    e007:	5d                   	pop    %rbp
    e008:	c3                   	ret    
    e009:	90                   	nop

000000000000e00a <_ZNK3MPI4Comm14Reduce_scatterEPKvPvPiRKNS_8DatatypeERKNS_2OpE>:
inline void
MPI::Comm::Reduce_scatter(const void *sendbuf, void *recvbuf,
                          int recvcounts[],
                          const MPI::Datatype & datatype,
                          const MPI::Op& op) const
{
    e00a:	f3 0f 1e fa          	endbr64 
    e00e:	55                   	push   %rbp
    e00f:	48 89 e5             	mov    %rsp,%rbp
    e012:	41 54                	push   %r12
    e014:	53                   	push   %rbx
    e015:	48 83 ec 30          	sub    $0x30,%rsp
    e019:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e01d:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e021:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    e025:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    e029:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    e02d:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    (void)MPI_Reduce_scatter(const_cast<void *>(sendbuf), recvbuf, recvcounts,
                             datatype, op, mpi_comm);
    e031:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e035:	4c 8b 60 08          	mov    0x8(%rax),%r12
    (void)MPI_Reduce_scatter(const_cast<void *>(sendbuf), recvbuf, recvcounts,
    e039:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    e03d:	48 89 c7             	mov    %rax,%rdi
    e040:	e8 af 1a 00 00       	call   faf4 <_ZNK3MPI2OpcvP9ompi_op_tEv>
    e045:	48 89 c3             	mov    %rax,%rbx
    e048:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e04c:	48 89 c7             	mov    %rax,%rdi
    e04f:	e8 2a d8 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    e054:	48 89 c1             	mov    %rax,%rcx
    e057:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    e05b:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    e05f:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e063:	4d 89 e1             	mov    %r12,%r9
    e066:	49 89 d8             	mov    %rbx,%r8
    e069:	48 89 c7             	mov    %rax,%rdi
    e06c:	e8 bf c4 ff ff       	call   a530 <MPI_Reduce_scatter@plt>
}
    e071:	90                   	nop
    e072:	48 83 c4 30          	add    $0x30,%rsp
    e076:	5b                   	pop    %rbx
    e077:	41 5c                	pop    %r12
    e079:	5d                   	pop    %rbp
    e07a:	c3                   	ret    
    e07b:	90                   	nop

000000000000e07c <_ZN3MPI4Comm10DisconnectEv>:
// Process Creation and Managemnt
//

inline void
MPI::Comm::Disconnect()
{
    e07c:	f3 0f 1e fa          	endbr64 
    e080:	55                   	push   %rbp
    e081:	48 89 e5             	mov    %rsp,%rbp
    e084:	48 83 ec 10          	sub    $0x10,%rsp
    e088:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
  (void) MPI_Comm_disconnect(&mpi_comm);
    e08c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e090:	48 83 c0 08          	add    $0x8,%rax
    e094:	48 89 c7             	mov    %rax,%rdi
    e097:	e8 44 be ff ff       	call   9ee0 <MPI_Comm_disconnect@plt>
}
    e09c:	90                   	nop
    e09d:	c9                   	leave  
    e09e:	c3                   	ret    
    e09f:	90                   	nop

000000000000e0a0 <_ZNK3MPI4Comm8Get_nameEPcRi>:
// External Interfaces
//

inline void
MPI::Comm::Get_name(char* comm_name, int& resultlen) const
{
    e0a0:	f3 0f 1e fa          	endbr64 
    e0a4:	55                   	push   %rbp
    e0a5:	48 89 e5             	mov    %rsp,%rbp
    e0a8:	48 83 ec 20          	sub    $0x20,%rsp
    e0ac:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    e0b0:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    e0b4:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
  (void) MPI_Comm_get_name(mpi_comm, comm_name, &resultlen);
    e0b8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e0bc:	48 8b 40 08          	mov    0x8(%rax),%rax
    e0c0:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    e0c4:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    e0c8:	48 89 ce             	mov    %rcx,%rsi
    e0cb:	48 89 c7             	mov    %rax,%rdi
    e0ce:	e8 6d bd ff ff       	call   9e40 <MPI_Comm_get_name@plt>
}
    e0d3:	90                   	nop
    e0d4:	c9                   	leave  
    e0d5:	c3                   	ret    

000000000000e0d6 <_ZN3MPI4Comm8Set_nameEPKc>:

inline void
MPI::Comm::Set_name(const char* comm_name)
{
    e0d6:	f3 0f 1e fa          	endbr64 
    e0da:	55                   	push   %rbp
    e0db:	48 89 e5             	mov    %rsp,%rbp
    e0de:	48 83 ec 10          	sub    $0x10,%rsp
    e0e2:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    e0e6:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
  (void) MPI_Comm_set_name(mpi_comm, const_cast<char *>(comm_name));
    e0ea:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e0ee:	48 8b 40 08          	mov    0x8(%rax),%rax
    e0f2:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    e0f6:	48 89 d6             	mov    %rdx,%rsi
    e0f9:	48 89 c7             	mov    %rax,%rdi
    e0fc:	e8 4f c8 ff ff       	call   a950 <MPI_Comm_set_name@plt>
}
    e101:	90                   	nop
    e102:	c9                   	leave  
    e103:	c3                   	ret    

000000000000e104 <_ZNK3MPI4Comm12Get_topologyEv>:
//Process Topologies
//

inline int
MPI::Comm::Get_topology() const
{
    e104:	f3 0f 1e fa          	endbr64 
    e108:	55                   	push   %rbp
    e109:	48 89 e5             	mov    %rsp,%rbp
    e10c:	48 83 ec 20          	sub    $0x20,%rsp
    e110:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e114:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e11b:	00 00 
    e11d:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e121:	31 c0                	xor    %eax,%eax
  int status;
  (void)MPI_Topo_test(mpi_comm, &status);
    e123:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e127:	48 8b 40 08          	mov    0x8(%rax),%rax
    e12b:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    e12f:	48 89 d6             	mov    %rdx,%rsi
    e132:	48 89 c7             	mov    %rax,%rdi
    e135:	e8 06 c6 ff ff       	call   a740 <MPI_Topo_test@plt>
  return status;
    e13a:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    e13d:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    e141:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    e148:	00 00 
    e14a:	74 05                	je     e151 <_ZNK3MPI4Comm12Get_topologyEv+0x4d>
    e14c:	e8 8f c3 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    e151:	c9                   	leave  
    e152:	c3                   	ret    
    e153:	90                   	nop

000000000000e154 <_ZN3MPI4Comm5AbortEi>:
// Environmental Inquiry
//

inline void
MPI::Comm::Abort(int errorcode)
{
    e154:	f3 0f 1e fa          	endbr64 
    e158:	55                   	push   %rbp
    e159:	48 89 e5             	mov    %rsp,%rbp
    e15c:	48 83 ec 10          	sub    $0x10,%rsp
    e160:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    e164:	89 75 f4             	mov    %esi,-0xc(%rbp)
  (void)MPI_Abort(mpi_comm, errorcode);
    e167:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e16b:	48 8b 40 08          	mov    0x8(%rax),%rax
    e16f:	8b 55 f4             	mov    -0xc(%rbp),%edx
    e172:	89 d6                	mov    %edx,%esi
    e174:	48 89 c7             	mov    %rax,%rdi
    e177:	e8 a4 c2 ff ff       	call   a420 <MPI_Abort@plt>
}
    e17c:	90                   	nop
    e17d:	c9                   	leave  
    e17e:	c3                   	ret    
    e17f:	90                   	nop

000000000000e180 <_ZNK3MPI4Comm14Get_errhandlerEv>:
//  going to be deprecated and replaced in MPI-2.
//

inline MPI::Errhandler
MPI::Comm::Get_errhandler() const
{
    e180:	f3 0f 1e fa          	endbr64 
    e184:	55                   	push   %rbp
    e185:	48 89 e5             	mov    %rsp,%rbp
    e188:	48 83 ec 20          	sub    $0x20,%rsp
    e18c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e190:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e194:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e19b:	00 00 
    e19d:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e1a1:	31 c0                	xor    %eax,%eax
    MPI_Errhandler errhandler;
    MPI_Comm_get_errhandler(mpi_comm, &errhandler);
    e1a3:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e1a7:	48 8b 40 08          	mov    0x8(%rax),%rax
    e1ab:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    e1af:	48 89 d6             	mov    %rdx,%rsi
    e1b2:	48 89 c7             	mov    %rax,%rdi
    e1b5:	e8 46 c4 ff ff       	call   a600 <MPI_Comm_get_errhandler@plt>
    return errhandler;
    e1ba:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    e1be:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e1c2:	48 89 d6             	mov    %rdx,%rsi
    e1c5:	48 89 c7             	mov    %rax,%rdi
    e1c8:	e8 b7 da ff ff       	call   bc84 <_ZN3MPI10ErrhandlerC1EP17ompi_errhandler_t>
}
    e1cd:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e1d1:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e1d8:	00 00 
    e1da:	74 05                	je     e1e1 <_ZNK3MPI4Comm14Get_errhandlerEv+0x61>
    e1dc:	e8 ff c2 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    e1e1:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e1e5:	c9                   	leave  
    e1e6:	c3                   	ret    
    e1e7:	90                   	nop

000000000000e1e8 <_ZN3MPI4Comm14Set_errhandlerERKNS_10ErrhandlerE>:

inline void
MPI::Comm::Set_errhandler(const MPI::Errhandler& errhandler)
{
    e1e8:	f3 0f 1e fa          	endbr64 
    e1ec:	55                   	push   %rbp
    e1ed:	48 89 e5             	mov    %rsp,%rbp
    e1f0:	48 83 ec 10          	sub    $0x10,%rsp
    e1f4:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    e1f8:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    (void)MPI_Comm_set_errhandler(mpi_comm, errhandler);
    e1fc:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    e200:	48 89 c7             	mov    %rax,%rdi
    e203:	e8 aa da ff ff       	call   bcb2 <_ZNK3MPI10ErrhandlercvP17ompi_errhandler_tEv>
    e208:	48 89 c2             	mov    %rax,%rdx
    e20b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e20f:	48 8b 40 08          	mov    0x8(%rax),%rax
    e213:	48 89 d6             	mov    %rdx,%rsi
    e216:	48 89 c7             	mov    %rax,%rdi
    e219:	e8 d2 c5 ff ff       	call   a7f0 <MPI_Comm_set_errhandler@plt>
}
    e21e:	90                   	nop
    e21f:	c9                   	leave  
    e220:	c3                   	ret    
    e221:	90                   	nop

000000000000e222 <_ZNK3MPI4Comm8Set_attrEiPKv>:
    (void) MPI_Comm_free_keyval(&comm_keyval);
}

inline void
MPI::Comm::Set_attr(int comm_keyval, const void* attribute_val) const
{
    e222:	f3 0f 1e fa          	endbr64 
    e226:	55                   	push   %rbp
    e227:	48 89 e5             	mov    %rsp,%rbp
    e22a:	48 83 ec 20          	sub    $0x20,%rsp
    e22e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    e232:	89 75 f4             	mov    %esi,-0xc(%rbp)
    e235:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    (void)MPI_Comm_set_attr(mpi_comm, comm_keyval, const_cast<void*>(attribute_val));
    e239:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e23d:	48 8b 40 08          	mov    0x8(%rax),%rax
    e241:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    e245:	8b 4d f4             	mov    -0xc(%rbp),%ecx
    e248:	89 ce                	mov    %ecx,%esi
    e24a:	48 89 c7             	mov    %rax,%rdi
    e24d:	e8 ee bd ff ff       	call   a040 <MPI_Comm_set_attr@plt>
}
    e252:	90                   	nop
    e253:	c9                   	leave  
    e254:	c3                   	ret    
    e255:	90                   	nop

000000000000e256 <_ZNK3MPI4Comm8Get_attrEiPv>:

inline bool
MPI::Comm::Get_attr(int comm_keyval, void* attribute_val) const
{
    e256:	f3 0f 1e fa          	endbr64 
    e25a:	55                   	push   %rbp
    e25b:	48 89 e5             	mov    %rsp,%rbp
    e25e:	48 83 ec 30          	sub    $0x30,%rsp
    e262:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e266:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    e269:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    e26d:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e274:	00 00 
    e276:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e27a:	31 c0                	xor    %eax,%eax
  int flag;
  (void)MPI_Comm_get_attr(mpi_comm, comm_keyval, attribute_val, &flag);
    e27c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e280:	48 8b 40 08          	mov    0x8(%rax),%rax
    e284:	48 8d 4d f4          	lea    -0xc(%rbp),%rcx
    e288:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
    e28c:	8b 75 e4             	mov    -0x1c(%rbp),%esi
    e28f:	48 89 c7             	mov    %rax,%rdi
    e292:	e8 49 c4 ff ff       	call   a6e0 <MPI_Comm_get_attr@plt>
  return OPAL_INT_TO_BOOL(flag);
    e297:	8b 45 f4             	mov    -0xc(%rbp),%eax
    e29a:	85 c0                	test   %eax,%eax
    e29c:	0f 95 c0             	setne  %al
}
    e29f:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    e2a3:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    e2aa:	00 00 
    e2ac:	74 05                	je     e2b3 <_ZNK3MPI4Comm8Get_attrEiPv+0x5d>
    e2ae:	e8 2d c2 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    e2b3:	c9                   	leave  
    e2b4:	c3                   	ret    
    e2b5:	90                   	nop

000000000000e2b6 <_ZN3MPI4Comm11Delete_attrEi>:

inline void
MPI::Comm::Delete_attr(int comm_keyval)
{
    e2b6:	f3 0f 1e fa          	endbr64 
    e2ba:	55                   	push   %rbp
    e2bb:	48 89 e5             	mov    %rsp,%rbp
    e2be:	48 83 ec 10          	sub    $0x10,%rsp
    e2c2:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    e2c6:	89 75 f4             	mov    %esi,-0xc(%rbp)
  (void)MPI_Comm_delete_attr(mpi_comm, comm_keyval);
    e2c9:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e2cd:	48 8b 40 08          	mov    0x8(%rax),%rax
    e2d1:	8b 55 f4             	mov    -0xc(%rbp),%edx
    e2d4:	89 d6                	mov    %edx,%esi
    e2d6:	48 89 c7             	mov    %rax,%rdi
    e2d9:	e8 b2 c1 ff ff       	call   a490 <MPI_Comm_delete_attr@plt>
}
    e2de:	90                   	nop
    e2df:	c9                   	leave  
    e2e0:	c3                   	ret    
    e2e1:	90                   	nop

000000000000e2e2 <_ZN3MPI9IntracommC1EP19ompi_communicator_t>:
//
// $HEADER$
//

inline
MPI::Intracomm::Intracomm(MPI_Comm data) {
    e2e2:	f3 0f 1e fa          	endbr64 
    e2e6:	55                   	push   %rbp
    e2e7:	48 89 e5             	mov    %rsp,%rbp
    e2ea:	53                   	push   %rbx
    e2eb:	48 83 ec 28          	sub    $0x28,%rsp
    e2ef:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    e2f3:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    e2f7:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e2fe:	00 00 
    e300:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    e304:	31 c0                	xor    %eax,%eax
    e306:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e30a:	48 89 c7             	mov    %rax,%rdi
    e30d:	e8 9e bd ff ff       	call   a0b0 <_ZN3MPI4CommC2Ev@plt>
    e312:	48 8d 15 3f db 00 00 	lea    0xdb3f(%rip),%rdx        # 1be58 <_ZTVN3MPI9IntracommE+0x10>
    e319:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e31d:	48 89 10             	mov    %rdx,(%rax)
  int flag = 0;
    e320:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
  if (MPI::Is_initialized() && (data != MPI_COMM_NULL)) {
    e327:	e8 a0 e5 ff ff       	call   c8cc <_ZN3MPI14Is_initializedEv>
    e32c:	84 c0                	test   %al,%al
    e32e:	74 14                	je     e344 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0x62>
    e330:	48 8d 05 e9 fa 00 00 	lea    0xfae9(%rip),%rax        # 1de20 <ompi_mpi_comm_null@@Base>
    e337:	48 39 45 d0          	cmp    %rax,-0x30(%rbp)
    e33b:	74 07                	je     e344 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0x62>
    e33d:	b8 01 00 00 00       	mov    $0x1,%eax
    e342:	eb 05                	jmp    e349 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0x67>
    e344:	b8 00 00 00 00       	mov    $0x0,%eax
    e349:	84 c0                	test   %al,%al
    e34b:	74 39                	je     e386 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0xa4>
    (void)MPI_Comm_test_inter(data, &flag);
    e34d:	48 8d 55 e4          	lea    -0x1c(%rbp),%rdx
    e351:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    e355:	48 89 d6             	mov    %rdx,%rsi
    e358:	48 89 c7             	mov    %rax,%rdi
    e35b:	e8 c0 c7 ff ff       	call   ab20 <MPI_Comm_test_inter@plt>
    if (flag) {
    e360:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    e363:	85 c0                	test   %eax,%eax
    e365:	74 11                	je     e378 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0x96>
      mpi_comm = MPI_COMM_NULL;
    e367:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e36b:	48 8d 15 ae fa 00 00 	lea    0xfaae(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    e372:	48 89 50 08          	mov    %rdx,0x8(%rax)
    }
  }
  else {
    mpi_comm = data;
  }
}
    e376:	eb 3a                	jmp    e3b2 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0xd0>
      mpi_comm = data;
    e378:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e37c:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    e380:	48 89 50 08          	mov    %rdx,0x8(%rax)
}
    e384:	eb 2c                	jmp    e3b2 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0xd0>
    mpi_comm = data;
    e386:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e38a:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    e38e:	48 89 50 08          	mov    %rdx,0x8(%rax)
}
    e392:	eb 1e                	jmp    e3b2 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0xd0>
    e394:	f3 0f 1e fa          	endbr64 
    e398:	48 89 c3             	mov    %rax,%rbx
    e39b:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e39f:	48 89 c7             	mov    %rax,%rdi
    e3a2:	e8 21 d9 ff ff       	call   bcc8 <_ZN3MPI4CommD1Ev>
    e3a7:	48 89 d8             	mov    %rbx,%rax
    e3aa:	48 89 c7             	mov    %rax,%rdi
    e3ad:	e8 fe c5 ff ff       	call   a9b0 <_Unwind_Resume@plt>
    e3b2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e3b6:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e3bd:	00 00 
    e3bf:	74 05                	je     e3c6 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0xe4>
    e3c1:	e8 1a c1 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    e3c6:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    e3ca:	c9                   	leave  
    e3cb:	c3                   	ret    

000000000000e3cc <_ZNK3MPI9Intracomm4ScanEPKvPviRKNS_8DatatypeERKNS_2OpE>:

inline void
MPI::Intracomm::Scan(const void *sendbuf, void *recvbuf, int count,
     const MPI::Datatype & datatype, const MPI::Op& op) const
{
    e3cc:	f3 0f 1e fa          	endbr64 
    e3d0:	55                   	push   %rbp
    e3d1:	48 89 e5             	mov    %rsp,%rbp
    e3d4:	41 54                	push   %r12
    e3d6:	53                   	push   %rbx
    e3d7:	48 83 ec 30          	sub    $0x30,%rsp
    e3db:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e3df:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e3e3:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    e3e7:	89 4d d4             	mov    %ecx,-0x2c(%rbp)
    e3ea:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    e3ee:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
  (void)MPI_Scan(const_cast<void *>(sendbuf), recvbuf, count, datatype, op, mpi_comm);
    e3f2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e3f6:	4c 8b 60 08          	mov    0x8(%rax),%r12
    e3fa:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    e3fe:	48 89 c7             	mov    %rax,%rdi
    e401:	e8 ee 16 00 00       	call   faf4 <_ZNK3MPI2OpcvP9ompi_op_tEv>
    e406:	48 89 c3             	mov    %rax,%rbx
    e409:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e40d:	48 89 c7             	mov    %rax,%rdi
    e410:	e8 69 d4 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    e415:	48 89 c1             	mov    %rax,%rcx
    e418:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    e41b:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    e41f:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e423:	4d 89 e1             	mov    %r12,%r9
    e426:	49 89 d8             	mov    %rbx,%r8
    e429:	48 89 c7             	mov    %rax,%rdi
    e42c:	e8 9f c0 ff ff       	call   a4d0 <MPI_Scan@plt>
}
    e431:	90                   	nop
    e432:	48 83 c4 30          	add    $0x30,%rsp
    e436:	5b                   	pop    %rbx
    e437:	41 5c                	pop    %r12
    e439:	5d                   	pop    %rbp
    e43a:	c3                   	ret    
    e43b:	90                   	nop

000000000000e43c <_ZNK3MPI9Intracomm6ExscanEPKvPviRKNS_8DatatypeERKNS_2OpE>:

inline void
MPI::Intracomm::Exscan(const void *sendbuf, void *recvbuf, int count,
			      const MPI::Datatype & datatype,
			      const MPI::Op& op) const
{
    e43c:	f3 0f 1e fa          	endbr64 
    e440:	55                   	push   %rbp
    e441:	48 89 e5             	mov    %rsp,%rbp
    e444:	41 54                	push   %r12
    e446:	53                   	push   %rbx
    e447:	48 83 ec 30          	sub    $0x30,%rsp
    e44b:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e44f:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e453:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    e457:	89 4d d4             	mov    %ecx,-0x2c(%rbp)
    e45a:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    e45e:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
  (void)MPI_Exscan(const_cast<void *>(sendbuf), recvbuf, count, datatype, op, mpi_comm);
    e462:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e466:	4c 8b 60 08          	mov    0x8(%rax),%r12
    e46a:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    e46e:	48 89 c7             	mov    %rax,%rdi
    e471:	e8 7e 16 00 00       	call   faf4 <_ZNK3MPI2OpcvP9ompi_op_tEv>
    e476:	48 89 c3             	mov    %rax,%rbx
    e479:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e47d:	48 89 c7             	mov    %rax,%rdi
    e480:	e8 f9 d3 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    e485:	48 89 c1             	mov    %rax,%rcx
    e488:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    e48b:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    e48f:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e493:	4d 89 e1             	mov    %r12,%r9
    e496:	49 89 d8             	mov    %rbx,%r8
    e499:	48 89 c7             	mov    %rax,%rdi
    e49c:	e8 ff be ff ff       	call   a3a0 <MPI_Exscan@plt>
}
    e4a1:	90                   	nop
    e4a2:	48 83 c4 30          	add    $0x30,%rsp
    e4a6:	5b                   	pop    %rbx
    e4a7:	41 5c                	pop    %r12
    e4a9:	5d                   	pop    %rbp
    e4aa:	c3                   	ret    
    e4ab:	90                   	nop

000000000000e4ac <_ZNK3MPI9Intracomm5CloneEv>:
  return newcomm;
}

inline MPI::Intracomm&
MPI::Intracomm::Clone() const
{
    e4ac:	f3 0f 1e fa          	endbr64 
    e4b0:	55                   	push   %rbp
    e4b1:	48 89 e5             	mov    %rsp,%rbp
    e4b4:	41 54                	push   %r12
    e4b6:	53                   	push   %rbx
    e4b7:	48 83 ec 30          	sub    $0x30,%rsp
    e4bb:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
    e4bf:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e4c6:	00 00 
    e4c8:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    e4cc:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void)MPI_Comm_dup(mpi_comm, &newcomm);
    e4ce:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e4d2:	48 8b 40 08          	mov    0x8(%rax),%rax
    e4d6:	48 8d 55 d8          	lea    -0x28(%rbp),%rdx
    e4da:	48 89 d6             	mov    %rdx,%rsi
    e4dd:	48 89 c7             	mov    %rax,%rdi
    e4e0:	e8 6b c1 ff ff       	call   a650 <MPI_Comm_dup@plt>
  MPI::Intracomm* dup = new MPI::Intracomm(newcomm);
    e4e5:	bf 10 00 00 00       	mov    $0x10,%edi
    e4ea:	e8 f1 be ff ff       	call   a3e0 <_Znwm@plt>
    e4ef:	48 89 c3             	mov    %rax,%rbx
    e4f2:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e4f6:	48 89 c6             	mov    %rax,%rsi
    e4f9:	48 89 df             	mov    %rbx,%rdi
    e4fc:	e8 e1 fd ff ff       	call   e2e2 <_ZN3MPI9IntracommC1EP19ompi_communicator_t>
    e501:	48 89 5d e0          	mov    %rbx,-0x20(%rbp)
  return *dup;
    e505:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
}
    e509:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    e50d:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    e514:	00 00 
    e516:	74 26                	je     e53e <_ZNK3MPI9Intracomm5CloneEv+0x92>
    e518:	eb 1f                	jmp    e539 <_ZNK3MPI9Intracomm5CloneEv+0x8d>
    e51a:	f3 0f 1e fa          	endbr64 
  MPI::Intracomm* dup = new MPI::Intracomm(newcomm);
    e51e:	49 89 c4             	mov    %rax,%r12
    e521:	be 10 00 00 00       	mov    $0x10,%esi
    e526:	48 89 df             	mov    %rbx,%rdi
    e529:	e8 d2 be ff ff       	call   a400 <_ZdlPvm@plt>
    e52e:	4c 89 e0             	mov    %r12,%rax
    e531:	48 89 c7             	mov    %rax,%rdi
    e534:	e8 77 c4 ff ff       	call   a9b0 <_Unwind_Resume@plt>
}
    e539:	e8 a2 bf ff ff       	call   a4e0 <__stack_chk_fail@plt>
    e53e:	48 83 c4 30          	add    $0x30,%rsp
    e542:	5b                   	pop    %rbx
    e543:	41 5c                	pop    %r12
    e545:	5d                   	pop    %rbp
    e546:	c3                   	ret    
    e547:	90                   	nop

000000000000e548 <_ZNK3MPI9Intracomm6CreateERKNS_5GroupE>:

inline MPI::Intracomm
MPI::Intracomm::Create(const MPI::Group& group) const
{
    e548:	f3 0f 1e fa          	endbr64 
    e54c:	55                   	push   %rbp
    e54d:	48 89 e5             	mov    %rsp,%rbp
    e550:	48 83 ec 30          	sub    $0x30,%rsp
    e554:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e558:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e55c:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    e560:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e567:	00 00 
    e569:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e56d:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void)MPI_Comm_create(mpi_comm, group, &newcomm);
    e56f:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e573:	48 89 c7             	mov    %rax,%rdi
    e576:	e8 8d d5 ff ff       	call   bb08 <_ZNK3MPI5GroupcvP12ompi_group_tEv>
    e57b:	48 89 c1             	mov    %rax,%rcx
    e57e:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e582:	48 8b 40 08          	mov    0x8(%rax),%rax
    e586:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    e58a:	48 89 ce             	mov    %rcx,%rsi
    e58d:	48 89 c7             	mov    %rax,%rdi
    e590:	e8 7b bd ff ff       	call   a310 <MPI_Comm_create@plt>
  return newcomm;
    e595:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    e599:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e59d:	48 89 d6             	mov    %rdx,%rsi
    e5a0:	48 89 c7             	mov    %rax,%rdi
    e5a3:	e8 3a fd ff ff       	call   e2e2 <_ZN3MPI9IntracommC1EP19ompi_communicator_t>
}
    e5a8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e5ac:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e5b3:	00 00 
    e5b5:	74 05                	je     e5bc <_ZNK3MPI9Intracomm6CreateERKNS_5GroupE+0x74>
    e5b7:	e8 24 bf ff ff       	call   a4e0 <__stack_chk_fail@plt>
    e5bc:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e5c0:	c9                   	leave  
    e5c1:	c3                   	ret    

000000000000e5c2 <_ZNK3MPI9Intracomm5SplitEii>:

inline MPI::Intracomm
MPI::Intracomm::Split(int color, int key) const
{
    e5c2:	f3 0f 1e fa          	endbr64 
    e5c6:	55                   	push   %rbp
    e5c7:	48 89 e5             	mov    %rsp,%rbp
    e5ca:	48 83 ec 30          	sub    $0x30,%rsp
    e5ce:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e5d2:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e5d6:	89 55 dc             	mov    %edx,-0x24(%rbp)
    e5d9:	89 4d d8             	mov    %ecx,-0x28(%rbp)
    e5dc:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e5e3:	00 00 
    e5e5:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e5e9:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void)MPI_Comm_split(mpi_comm, color, key, &newcomm);
    e5eb:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e5ef:	48 8b 40 08          	mov    0x8(%rax),%rax
    e5f3:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    e5f7:	8b 55 d8             	mov    -0x28(%rbp),%edx
    e5fa:	8b 75 dc             	mov    -0x24(%rbp),%esi
    e5fd:	48 89 c7             	mov    %rax,%rdi
    e600:	e8 8b c2 ff ff       	call   a890 <MPI_Comm_split@plt>
  return newcomm;
    e605:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    e609:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e60d:	48 89 d6             	mov    %rdx,%rsi
    e610:	48 89 c7             	mov    %rax,%rdi
    e613:	e8 ca fc ff ff       	call   e2e2 <_ZN3MPI9IntracommC1EP19ompi_communicator_t>
}
    e618:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e61c:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e623:	00 00 
    e625:	74 05                	je     e62c <_ZNK3MPI9Intracomm5SplitEii+0x6a>
    e627:	e8 b4 be ff ff       	call   a4e0 <__stack_chk_fail@plt>
    e62c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e630:	c9                   	leave  
    e631:	c3                   	ret    

000000000000e632 <_ZNK3MPI9Intracomm16Create_intercommEiRKNS_4CommEii>:

inline MPI::Intercomm
MPI::Intracomm::Create_intercomm(int local_leader,
					const MPI::Comm& peer_comm,
					int remote_leader, int tag) const
{
    e632:	f3 0f 1e fa          	endbr64 
    e636:	55                   	push   %rbp
    e637:	48 89 e5             	mov    %rsp,%rbp
    e63a:	48 83 ec 40          	sub    $0x40,%rsp
    e63e:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e642:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e646:	89 55 dc             	mov    %edx,-0x24(%rbp)
    e649:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    e64d:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    e651:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
    e655:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e65c:	00 00 
    e65e:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e662:	31 c0                	xor    %eax,%eax
  MPI_Comm newintercomm;
  (void)MPI_Intercomm_create(mpi_comm, local_leader, peer_comm,
    e664:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    e668:	48 89 c7             	mov    %rax,%rdi
    e66b:	e8 2a d5 ff ff       	call   bb9a <_ZNK3MPI9Comm_NullcvP19ompi_communicator_tEv>
    e670:	48 89 c7             	mov    %rax,%rdi
    e673:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e677:	48 8b 40 08          	mov    0x8(%rax),%rax
    e67b:	4c 8d 45 f0          	lea    -0x10(%rbp),%r8
    e67f:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    e682:	8b 55 d8             	mov    -0x28(%rbp),%edx
    e685:	8b 75 dc             	mov    -0x24(%rbp),%esi
    e688:	4d 89 c1             	mov    %r8,%r9
    e68b:	41 89 c8             	mov    %ecx,%r8d
    e68e:	89 d1                	mov    %edx,%ecx
    e690:	48 89 fa             	mov    %rdi,%rdx
    e693:	48 89 c7             	mov    %rax,%rdi
    e696:	e8 95 bd ff ff       	call   a430 <MPI_Intercomm_create@plt>
			     remote_leader, tag, &newintercomm);
  return newintercomm;
    e69b:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    e69f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e6a3:	48 89 d6             	mov    %rdx,%rsi
    e6a6:	48 89 c7             	mov    %rax,%rdi
    e6a9:	e8 04 d7 ff ff       	call   bdb2 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
}
    e6ae:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e6b2:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e6b9:	00 00 
    e6bb:	74 05                	je     e6c2 <_ZNK3MPI9Intracomm16Create_intercommEiRKNS_4CommEii+0x90>
    e6bd:	e8 1e be ff ff       	call   a4e0 <__stack_chk_fail@plt>
    e6c2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e6c6:	c9                   	leave  
    e6c7:	c3                   	ret    

000000000000e6c8 <_ZNK3MPI9Intracomm11Create_cartEiPKiPKbb>:

inline MPI::Cartcomm
MPI::Intracomm::Create_cart(int ndims, const int dims[],
				   const bool periods[], bool reorder) const
{
    e6c8:	f3 0f 1e fa          	endbr64 
    e6cc:	55                   	push   %rbp
    e6cd:	48 89 e5             	mov    %rsp,%rbp
    e6d0:	48 83 ec 50          	sub    $0x50,%rsp
    e6d4:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    e6d8:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    e6dc:	89 55 cc             	mov    %edx,-0x34(%rbp)
    e6df:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    e6e3:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    e6e7:	44 89 c8             	mov    %r9d,%eax
    e6ea:	88 45 c8             	mov    %al,-0x38(%rbp)
    e6ed:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e6f4:	00 00 
    e6f6:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e6fa:	31 c0                	xor    %eax,%eax
  int *int_periods = new int [ndims];
    e6fc:	8b 45 cc             	mov    -0x34(%rbp),%eax
    e6ff:	48 98                	cltq   
    e701:	48 ba fe ff ff ff ff 	movabs $0x1ffffffffffffffe,%rdx
    e708:	ff ff 1f 
    e70b:	48 39 d0             	cmp    %rdx,%rax
    e70e:	77 19                	ja     e729 <_ZNK3MPI9Intracomm11Create_cartEiPKiPKbb+0x61>
    e710:	48 c1 e0 02          	shl    $0x2,%rax
    e714:	48 89 c7             	mov    %rax,%rdi
    e717:	e8 a4 b6 ff ff       	call   9dc0 <_Znam@plt>
    e71c:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
  for (int i=0; i<ndims; i++)
    e720:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
    e727:	eb 33                	jmp    e75c <_ZNK3MPI9Intracomm11Create_cartEiPKiPKbb+0x94>
  int *int_periods = new int [ndims];
    e729:	e8 c2 bf ff ff       	call   a6f0 <__cxa_throw_bad_array_new_length@plt>
    int_periods[i] = (int) periods[i];
    e72e:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    e731:	48 63 d0             	movslq %eax,%rdx
    e734:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    e738:	48 01 d0             	add    %rdx,%rax
    e73b:	0f b6 00             	movzbl (%rax),%eax
    e73e:	8b 55 e4             	mov    -0x1c(%rbp),%edx
    e741:	48 63 d2             	movslq %edx,%rdx
    e744:	48 8d 0c 95 00 00 00 	lea    0x0(,%rdx,4),%rcx
    e74b:	00 
    e74c:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    e750:	48 01 ca             	add    %rcx,%rdx
    e753:	0f b6 c0             	movzbl %al,%eax
    e756:	89 02                	mov    %eax,(%rdx)
  for (int i=0; i<ndims; i++)
    e758:	83 45 e4 01          	addl   $0x1,-0x1c(%rbp)
    e75c:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    e75f:	3b 45 cc             	cmp    -0x34(%rbp),%eax
    e762:	7c ca                	jl     e72e <_ZNK3MPI9Intracomm11Create_cartEiPKiPKbb+0x66>

  MPI_Comm newcomm;
  (void)MPI_Cart_create(mpi_comm, ndims, const_cast<int *>(dims),
    e764:	0f b6 7d c8          	movzbl -0x38(%rbp),%edi
    e768:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    e76c:	48 8b 40 08          	mov    0x8(%rax),%rax
    e770:	4c 8d 45 e8          	lea    -0x18(%rbp),%r8
    e774:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    e778:	48 8b 55 c0          	mov    -0x40(%rbp),%rdx
    e77c:	8b 75 cc             	mov    -0x34(%rbp),%esi
    e77f:	4d 89 c1             	mov    %r8,%r9
    e782:	41 89 f8             	mov    %edi,%r8d
    e785:	48 89 c7             	mov    %rax,%rdi
    e788:	e8 a3 b6 ff ff       	call   9e30 <MPI_Cart_create@plt>
		      int_periods, (int)reorder, &newcomm);
  delete [] int_periods;
    e78d:	48 83 7d f0 00       	cmpq   $0x0,-0x10(%rbp)
    e792:	74 0c                	je     e7a0 <_ZNK3MPI9Intracomm11Create_cartEiPKiPKbb+0xd8>
    e794:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    e798:	48 89 c7             	mov    %rax,%rdi
    e79b:	e8 30 be ff ff       	call   a5d0 <_ZdaPv@plt>
  return newcomm;
    e7a0:	48 8d 55 e8          	lea    -0x18(%rbp),%rdx
    e7a4:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e7a8:	48 89 d6             	mov    %rdx,%rsi
    e7ab:	48 89 c7             	mov    %rax,%rdi
    e7ae:	e8 47 05 00 00       	call   ecfa <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t>
}
    e7b3:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e7b7:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e7be:	00 00 
    e7c0:	74 05                	je     e7c7 <_ZNK3MPI9Intracomm11Create_cartEiPKiPKbb+0xff>
    e7c2:	e8 19 bd ff ff       	call   a4e0 <__stack_chk_fail@plt>
    e7c7:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e7cb:	c9                   	leave  
    e7cc:	c3                   	ret    
    e7cd:	90                   	nop

000000000000e7ce <_ZNK3MPI9Intracomm12Create_graphEiPKiS2_b>:

inline MPI::Graphcomm
MPI::Intracomm::Create_graph(int nnodes, const int index[],
				    const int edges[], bool reorder) const
{
    e7ce:	f3 0f 1e fa          	endbr64 
    e7d2:	55                   	push   %rbp
    e7d3:	48 89 e5             	mov    %rsp,%rbp
    e7d6:	48 83 ec 40          	sub    $0x40,%rsp
    e7da:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e7de:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e7e2:	89 55 dc             	mov    %edx,-0x24(%rbp)
    e7e5:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    e7e9:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    e7ed:	44 89 c8             	mov    %r9d,%eax
    e7f0:	88 45 d8             	mov    %al,-0x28(%rbp)
    e7f3:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e7fa:	00 00 
    e7fc:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e800:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void)MPI_Graph_create(mpi_comm, nnodes, const_cast<int *>(index),
    e802:	0f b6 7d d8          	movzbl -0x28(%rbp),%edi
    e806:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e80a:	48 8b 40 08          	mov    0x8(%rax),%rax
    e80e:	4c 8d 45 f0          	lea    -0x10(%rbp),%r8
    e812:	48 8b 4d c8          	mov    -0x38(%rbp),%rcx
    e816:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    e81a:	8b 75 dc             	mov    -0x24(%rbp),%esi
    e81d:	4d 89 c1             	mov    %r8,%r9
    e820:	41 89 f8             	mov    %edi,%r8d
    e823:	48 89 c7             	mov    %rax,%rdi
    e826:	e8 e5 c2 ff ff       	call   ab10 <MPI_Graph_create@plt>
                         const_cast<int *>(edges), (int)reorder, &newcomm);
  return newcomm;
    e82b:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    e82f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e833:	48 89 d6             	mov    %rdx,%rsi
    e836:	48 89 c7             	mov    %rax,%rdi
    e839:	e8 54 0a 00 00       	call   f292 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t>
}
    e83e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e842:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e849:	00 00 
    e84b:	74 05                	je     e852 <_ZNK3MPI9Intracomm12Create_graphEiPKiS2_b+0x84>
    e84d:	e8 8e bc ff ff       	call   a4e0 <__stack_chk_fail@plt>
    e852:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e856:	c9                   	leave  
    e857:	c3                   	ret    

000000000000e858 <_ZNK3MPI9Intracomm6AcceptEPKcRKNS_4InfoEi>:

inline MPI::Intercomm
MPI::Intracomm::Accept(const char* port_name,
			      const MPI::Info& info,
			      int root) const
{
    e858:	f3 0f 1e fa          	endbr64 
    e85c:	55                   	push   %rbp
    e85d:	48 89 e5             	mov    %rsp,%rbp
    e860:	53                   	push   %rbx
    e861:	48 83 ec 48          	sub    $0x48,%rsp
    e865:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    e869:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    e86d:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    e871:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    e875:	44 89 45 bc          	mov    %r8d,-0x44(%rbp)
    e879:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e880:	00 00 
    e882:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    e886:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void) MPI_Comm_accept(const_cast<char *>(port_name), info, root, mpi_comm,
    e888:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    e88c:	48 8b 58 08          	mov    0x8(%rax),%rbx
    e890:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    e894:	48 89 c7             	mov    %rax,%rdi
    e897:	e8 ca d5 ff ff       	call   be66 <_ZNK3MPI4InfocvP11ompi_info_tEv>
    e89c:	48 89 c6             	mov    %rax,%rsi
    e89f:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    e8a3:	8b 55 bc             	mov    -0x44(%rbp),%edx
    e8a6:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e8aa:	49 89 c8             	mov    %rcx,%r8
    e8ad:	48 89 d9             	mov    %rbx,%rcx
    e8b0:	48 89 c7             	mov    %rax,%rdi
    e8b3:	e8 28 c1 ff ff       	call   a9e0 <MPI_Comm_accept@plt>
			 &newcomm);
  return newcomm;
    e8b8:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    e8bc:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e8c0:	48 89 d6             	mov    %rdx,%rsi
    e8c3:	48 89 c7             	mov    %rax,%rdi
    e8c6:	e8 e7 d4 ff ff       	call   bdb2 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
}
    e8cb:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e8cf:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e8d6:	00 00 
    e8d8:	74 05                	je     e8df <_ZNK3MPI9Intracomm6AcceptEPKcRKNS_4InfoEi+0x87>
    e8da:	e8 01 bc ff ff       	call   a4e0 <__stack_chk_fail@plt>
    e8df:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e8e3:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    e8e7:	c9                   	leave  
    e8e8:	c3                   	ret    
    e8e9:	90                   	nop

000000000000e8ea <_ZNK3MPI9Intracomm7ConnectEPKcRKNS_4InfoEi>:

inline MPI::Intercomm
MPI::Intracomm::Connect(const char* port_name,
			       const MPI::Info& info,
			       int root) const
{
    e8ea:	f3 0f 1e fa          	endbr64 
    e8ee:	55                   	push   %rbp
    e8ef:	48 89 e5             	mov    %rsp,%rbp
    e8f2:	53                   	push   %rbx
    e8f3:	48 83 ec 48          	sub    $0x48,%rsp
    e8f7:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    e8fb:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    e8ff:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    e903:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    e907:	44 89 45 bc          	mov    %r8d,-0x44(%rbp)
    e90b:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e912:	00 00 
    e914:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    e918:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void) MPI_Comm_connect(const_cast<char *>(port_name), info, root, mpi_comm,
    e91a:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    e91e:	48 8b 58 08          	mov    0x8(%rax),%rbx
    e922:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    e926:	48 89 c7             	mov    %rax,%rdi
    e929:	e8 38 d5 ff ff       	call   be66 <_ZNK3MPI4InfocvP11ompi_info_tEv>
    e92e:	48 89 c6             	mov    %rax,%rsi
    e931:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    e935:	8b 55 bc             	mov    -0x44(%rbp),%edx
    e938:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e93c:	49 89 c8             	mov    %rcx,%r8
    e93f:	48 89 d9             	mov    %rbx,%rcx
    e942:	48 89 c7             	mov    %rax,%rdi
    e945:	e8 f6 be ff ff       	call   a840 <MPI_Comm_connect@plt>
			  &newcomm);
  return newcomm;
    e94a:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    e94e:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e952:	48 89 d6             	mov    %rdx,%rsi
    e955:	48 89 c7             	mov    %rax,%rdi
    e958:	e8 55 d4 ff ff       	call   bdb2 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
}
    e95d:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e961:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e968:	00 00 
    e96a:	74 05                	je     e971 <_ZNK3MPI9Intracomm7ConnectEPKcRKNS_4InfoEi+0x87>
    e96c:	e8 6f bb ff ff       	call   a4e0 <__stack_chk_fail@plt>
    e971:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e975:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    e979:	c9                   	leave  
    e97a:	c3                   	ret    
    e97b:	90                   	nop

000000000000e97c <_ZNK3MPI9Intracomm5SpawnEPKcPS2_iRKNS_4InfoEi>:

inline MPI::Intercomm
MPI::Intracomm::Spawn(const char* command, const char* argv[],
			     int maxprocs, const MPI::Info& info,
			     int root) const
{
    e97c:	f3 0f 1e fa          	endbr64 
    e980:	55                   	push   %rbp
    e981:	48 89 e5             	mov    %rsp,%rbp
    e984:	53                   	push   %rbx
    e985:	48 83 ec 48          	sub    $0x48,%rsp
    e989:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    e98d:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    e991:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    e995:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    e999:	44 89 45 bc          	mov    %r8d,-0x44(%rbp)
    e99d:	4c 89 4d b0          	mov    %r9,-0x50(%rbp)
    e9a1:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e9a8:	00 00 
    e9aa:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    e9ae:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void) MPI_Comm_spawn(const_cast<char *>(command), const_cast<char **>(argv), maxprocs,
			info, root, mpi_comm, &newcomm,
    e9b0:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    e9b4:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void) MPI_Comm_spawn(const_cast<char *>(command), const_cast<char **>(argv), maxprocs,
    e9b8:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
    e9bc:	48 89 c7             	mov    %rax,%rdi
    e9bf:	e8 a2 d4 ff ff       	call   be66 <_ZNK3MPI4InfocvP11ompi_info_tEv>
    e9c4:	48 89 c7             	mov    %rax,%rdi
    e9c7:	8b 55 bc             	mov    -0x44(%rbp),%edx
    e9ca:	48 8b 75 c0          	mov    -0x40(%rbp),%rsi
    e9ce:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e9d2:	6a 00                	push   $0x0
    e9d4:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    e9d8:	51                   	push   %rcx
    e9d9:	49 89 d9             	mov    %rbx,%r9
    e9dc:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    e9e0:	48 89 f9             	mov    %rdi,%rcx
    e9e3:	48 89 c7             	mov    %rax,%rdi
    e9e6:	e8 d5 b7 ff ff       	call   a1c0 <MPI_Comm_spawn@plt>
    e9eb:	48 83 c4 10          	add    $0x10,%rsp
			(int *)MPI_ERRCODES_IGNORE);
  return newcomm;
    e9ef:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    e9f3:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e9f7:	48 89 d6             	mov    %rdx,%rsi
    e9fa:	48 89 c7             	mov    %rax,%rdi
    e9fd:	e8 b0 d3 ff ff       	call   bdb2 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
}
    ea02:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ea06:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    ea0d:	00 00 
    ea0f:	74 05                	je     ea16 <_ZNK3MPI9Intracomm5SpawnEPKcPS2_iRKNS_4InfoEi+0x9a>
    ea11:	e8 ca ba ff ff       	call   a4e0 <__stack_chk_fail@plt>
    ea16:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ea1a:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    ea1e:	c9                   	leave  
    ea1f:	c3                   	ret    

000000000000ea20 <_ZNK3MPI9Intracomm5SpawnEPKcPS2_iRKNS_4InfoEiPi>:

inline MPI::Intercomm
MPI::Intracomm::Spawn(const char* command, const char* argv[],
                             int maxprocs, const MPI::Info& info,
                             int root, int array_of_errcodes[]) const
{
    ea20:	f3 0f 1e fa          	endbr64 
    ea24:	55                   	push   %rbp
    ea25:	48 89 e5             	mov    %rsp,%rbp
    ea28:	53                   	push   %rbx
    ea29:	48 83 ec 58          	sub    $0x58,%rsp
    ea2d:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    ea31:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    ea35:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    ea39:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    ea3d:	44 89 45 bc          	mov    %r8d,-0x44(%rbp)
    ea41:	4c 89 4d b0          	mov    %r9,-0x50(%rbp)
    ea45:	48 8b 45 18          	mov    0x18(%rbp),%rax
    ea49:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
    ea4d:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ea54:	00 00 
    ea56:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    ea5a:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void) MPI_Comm_spawn(const_cast<char *>(command), const_cast<char **>(argv), maxprocs,
                        info, root, mpi_comm, &newcomm,
    ea5c:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ea60:	48 8b 58 08          	mov    0x8(%rax),%rbx
  (void) MPI_Comm_spawn(const_cast<char *>(command), const_cast<char **>(argv), maxprocs,
    ea64:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
    ea68:	48 89 c7             	mov    %rax,%rdi
    ea6b:	e8 f6 d3 ff ff       	call   be66 <_ZNK3MPI4InfocvP11ompi_info_tEv>
    ea70:	48 89 c7             	mov    %rax,%rdi
    ea73:	8b 55 bc             	mov    -0x44(%rbp),%edx
    ea76:	48 8b 75 c0          	mov    -0x40(%rbp),%rsi
    ea7a:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    ea7e:	ff 75 a8             	push   -0x58(%rbp)
    ea81:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    ea85:	51                   	push   %rcx
    ea86:	49 89 d9             	mov    %rbx,%r9
    ea89:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    ea8d:	48 89 f9             	mov    %rdi,%rcx
    ea90:	48 89 c7             	mov    %rax,%rdi
    ea93:	e8 28 b7 ff ff       	call   a1c0 <MPI_Comm_spawn@plt>
    ea98:	48 83 c4 10          	add    $0x10,%rsp
			array_of_errcodes);
  return newcomm;
    ea9c:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    eaa0:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    eaa4:	48 89 d6             	mov    %rdx,%rsi
    eaa7:	48 89 c7             	mov    %rax,%rdi
    eaaa:	e8 03 d3 ff ff       	call   bdb2 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
}
    eaaf:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    eab3:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    eaba:	00 00 
    eabc:	74 05                	je     eac3 <_ZNK3MPI9Intracomm5SpawnEPKcPS2_iRKNS_4InfoEiPi+0xa3>
    eabe:	e8 1d ba ff ff       	call   a4e0 <__stack_chk_fail@plt>
    eac3:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    eac7:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    eacb:	c9                   	leave  
    eacc:	c3                   	ret    
    eacd:	90                   	nop

000000000000eace <_ZN3MPI9Intracomm14Spawn_multipleEiPPKcPS3_PKiPKNS_4InfoEi>:
MPI::Intracomm::Spawn_multiple(int count,
				      const char* array_of_commands[],
				      const char** array_of_argv[],
				      const int array_of_maxprocs[],
				      const Info array_of_info[], int root)
{
    eace:	f3 0f 1e fa          	endbr64 
    ead2:	55                   	push   %rbp
    ead3:	48 89 e5             	mov    %rsp,%rbp
    ead6:	48 83 ec 60          	sub    $0x60,%rsp
    eada:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    eade:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    eae2:	89 55 cc             	mov    %edx,-0x34(%rbp)
    eae5:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    eae9:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    eaed:	4c 89 4d b0          	mov    %r9,-0x50(%rbp)
    eaf1:	48 8b 45 10          	mov    0x10(%rbp),%rax
    eaf5:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
    eaf9:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    eb00:	00 00 
    eb02:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    eb06:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  MPI_Info *const array_of_mpi_info =
      convert_info_to_mpi_info(count, array_of_info);
    eb08:	48 8b 55 a8          	mov    -0x58(%rbp),%rdx
    eb0c:	8b 45 cc             	mov    -0x34(%rbp),%eax
    eb0f:	48 89 d6             	mov    %rdx,%rsi
    eb12:	89 c7                	mov    %eax,%edi
    eb14:	e8 7f 00 00 00       	call   eb98 <_ZN3MPI9Intracomm24convert_info_to_mpi_infoEiPKNS_4InfoE>
    eb19:	48 89 45 f0          	mov    %rax,-0x10(%rbp)

  MPI_Comm_spawn_multiple(count, const_cast<char **>(array_of_commands),
    eb1d:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    eb21:	48 8b 78 08          	mov    0x8(%rax),%rdi
    eb25:	44 8b 4d 18          	mov    0x18(%rbp),%r9d
    eb29:	4c 8b 55 f0          	mov    -0x10(%rbp),%r10
    eb2d:	48 8b 4d b0          	mov    -0x50(%rbp),%rcx
    eb31:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
    eb35:	48 8b 75 c0          	mov    -0x40(%rbp),%rsi
    eb39:	8b 45 cc             	mov    -0x34(%rbp),%eax
    eb3c:	48 83 ec 08          	sub    $0x8,%rsp
    eb40:	6a 00                	push   $0x0
    eb42:	4c 8d 45 e8          	lea    -0x18(%rbp),%r8
    eb46:	41 50                	push   %r8
    eb48:	57                   	push   %rdi
    eb49:	4d 89 d0             	mov    %r10,%r8
    eb4c:	89 c7                	mov    %eax,%edi
    eb4e:	e8 6d bd ff ff       	call   a8c0 <MPI_Comm_spawn_multiple@plt>
    eb53:	48 83 c4 20          	add    $0x20,%rsp
			  const_cast<char ***>(array_of_argv),
                          const_cast<int *>(array_of_maxprocs),
			  array_of_mpi_info, root,
			  mpi_comm, &newcomm, (int *)MPI_ERRCODES_IGNORE);
  delete[] array_of_mpi_info;
    eb57:	48 83 7d f0 00       	cmpq   $0x0,-0x10(%rbp)
    eb5c:	74 0c                	je     eb6a <_ZN3MPI9Intracomm14Spawn_multipleEiPPKcPS3_PKiPKNS_4InfoEi+0x9c>
    eb5e:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    eb62:	48 89 c7             	mov    %rax,%rdi
    eb65:	e8 66 ba ff ff       	call   a5d0 <_ZdaPv@plt>
  return newcomm;
    eb6a:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    eb6e:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    eb72:	48 89 d6             	mov    %rdx,%rsi
    eb75:	48 89 c7             	mov    %rax,%rdi
    eb78:	e8 35 d2 ff ff       	call   bdb2 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
}
    eb7d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    eb81:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    eb88:	00 00 
    eb8a:	74 05                	je     eb91 <_ZN3MPI9Intracomm14Spawn_multipleEiPPKcPS3_PKiPKNS_4InfoEi+0xc3>
    eb8c:	e8 4f b9 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    eb91:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    eb95:	c9                   	leave  
    eb96:	c3                   	ret    
    eb97:	90                   	nop

000000000000eb98 <_ZN3MPI9Intracomm24convert_info_to_mpi_infoEiPKNS_4InfoE>:

inline MPI_Info *
MPI::Intracomm::convert_info_to_mpi_info(int p_nbr, const Info p_info_tbl[])
{
    eb98:	f3 0f 1e fa          	endbr64 
    eb9c:	55                   	push   %rbp
    eb9d:	48 89 e5             	mov    %rsp,%rbp
    eba0:	53                   	push   %rbx
    eba1:	48 83 ec 28          	sub    $0x28,%rsp
    eba5:	89 7d dc             	mov    %edi,-0x24(%rbp)
    eba8:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
   MPI_Info *const mpi_info_tbl = new MPI_Info [p_nbr];
    ebac:	8b 45 dc             	mov    -0x24(%rbp),%eax
    ebaf:	48 98                	cltq   
    ebb1:	48 ba ff ff ff ff ff 	movabs $0xfffffffffffffff,%rdx
    ebb8:	ff ff 0f 
    ebbb:	48 39 d0             	cmp    %rdx,%rax
    ebbe:	77 19                	ja     ebd9 <_ZN3MPI9Intracomm24convert_info_to_mpi_infoEiPKNS_4InfoE+0x41>
    ebc0:	48 c1 e0 03          	shl    $0x3,%rax
    ebc4:	48 89 c7             	mov    %rax,%rdi
    ebc7:	e8 f4 b1 ff ff       	call   9dc0 <_Znam@plt>
    ebcc:	48 89 45 e8          	mov    %rax,-0x18(%rbp)

   for (int i_tbl=0; i_tbl < p_nbr; i_tbl++) {
    ebd0:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
    ebd7:	eb 3d                	jmp    ec16 <_ZN3MPI9Intracomm24convert_info_to_mpi_infoEiPKNS_4InfoE+0x7e>
   MPI_Info *const mpi_info_tbl = new MPI_Info [p_nbr];
    ebd9:	e8 12 bb ff ff       	call   a6f0 <__cxa_throw_bad_array_new_length@plt>
       mpi_info_tbl[i_tbl] = p_info_tbl[i_tbl];
    ebde:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    ebe1:	48 98                	cltq   
    ebe3:	48 c1 e0 04          	shl    $0x4,%rax
    ebe7:	48 89 c2             	mov    %rax,%rdx
    ebea:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ebee:	48 01 d0             	add    %rdx,%rax
    ebf1:	8b 55 e4             	mov    -0x1c(%rbp),%edx
    ebf4:	48 63 d2             	movslq %edx,%rdx
    ebf7:	48 8d 0c d5 00 00 00 	lea    0x0(,%rdx,8),%rcx
    ebfe:	00 
    ebff:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    ec03:	48 8d 1c 11          	lea    (%rcx,%rdx,1),%rbx
    ec07:	48 89 c7             	mov    %rax,%rdi
    ec0a:	e8 57 d2 ff ff       	call   be66 <_ZNK3MPI4InfocvP11ompi_info_tEv>
    ec0f:	48 89 03             	mov    %rax,(%rbx)
   for (int i_tbl=0; i_tbl < p_nbr; i_tbl++) {
    ec12:	83 45 e4 01          	addl   $0x1,-0x1c(%rbp)
    ec16:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    ec19:	3b 45 dc             	cmp    -0x24(%rbp),%eax
    ec1c:	7c c0                	jl     ebde <_ZN3MPI9Intracomm24convert_info_to_mpi_infoEiPKNS_4InfoE+0x46>
   }

   return mpi_info_tbl;
    ec1e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
}
    ec22:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    ec26:	c9                   	leave  
    ec27:	c3                   	ret    

000000000000ec28 <_ZN3MPI9Intracomm14Spawn_multipleEiPPKcPS3_PKiPKNS_4InfoEiPi>:
                                      const char* array_of_commands[],
                                      const char** array_of_argv[],
                                      const int array_of_maxprocs[],
                                      const Info array_of_info[], int root,
				      int array_of_errcodes[])
{
    ec28:	f3 0f 1e fa          	endbr64 
    ec2c:	55                   	push   %rbp
    ec2d:	48 89 e5             	mov    %rsp,%rbp
    ec30:	48 83 ec 60          	sub    $0x60,%rsp
    ec34:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    ec38:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    ec3c:	89 55 cc             	mov    %edx,-0x34(%rbp)
    ec3f:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    ec43:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    ec47:	4c 89 4d b0          	mov    %r9,-0x50(%rbp)
    ec4b:	48 8b 45 10          	mov    0x10(%rbp),%rax
    ec4f:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
    ec53:	48 8b 45 20          	mov    0x20(%rbp),%rax
    ec57:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
    ec5b:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ec62:	00 00 
    ec64:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    ec68:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  MPI_Info *const array_of_mpi_info =
      convert_info_to_mpi_info(count, array_of_info);
    ec6a:	48 8b 55 a8          	mov    -0x58(%rbp),%rdx
    ec6e:	8b 45 cc             	mov    -0x34(%rbp),%eax
    ec71:	48 89 d6             	mov    %rdx,%rsi
    ec74:	89 c7                	mov    %eax,%edi
    ec76:	e8 1d ff ff ff       	call   eb98 <_ZN3MPI9Intracomm24convert_info_to_mpi_infoEiPKNS_4InfoE>
    ec7b:	48 89 45 f0          	mov    %rax,-0x10(%rbp)

  MPI_Comm_spawn_multiple(count, const_cast<char **>(array_of_commands),
    ec7f:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ec83:	48 8b 78 08          	mov    0x8(%rax),%rdi
    ec87:	44 8b 4d 18          	mov    0x18(%rbp),%r9d
    ec8b:	4c 8b 55 f0          	mov    -0x10(%rbp),%r10
    ec8f:	48 8b 4d b0          	mov    -0x50(%rbp),%rcx
    ec93:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
    ec97:	48 8b 75 c0          	mov    -0x40(%rbp),%rsi
    ec9b:	8b 45 cc             	mov    -0x34(%rbp),%eax
    ec9e:	48 83 ec 08          	sub    $0x8,%rsp
    eca2:	ff 75 a0             	push   -0x60(%rbp)
    eca5:	4c 8d 45 e8          	lea    -0x18(%rbp),%r8
    eca9:	41 50                	push   %r8
    ecab:	57                   	push   %rdi
    ecac:	4d 89 d0             	mov    %r10,%r8
    ecaf:	89 c7                	mov    %eax,%edi
    ecb1:	e8 0a bc ff ff       	call   a8c0 <MPI_Comm_spawn_multiple@plt>
    ecb6:	48 83 c4 20          	add    $0x20,%rsp
                          const_cast<char ***>(array_of_argv),
                          const_cast<int *>(array_of_maxprocs),
                          array_of_mpi_info, root,
                          mpi_comm, &newcomm, array_of_errcodes);
  delete[] array_of_mpi_info;
    ecba:	48 83 7d f0 00       	cmpq   $0x0,-0x10(%rbp)
    ecbf:	74 0c                	je     eccd <_ZN3MPI9Intracomm14Spawn_multipleEiPPKcPS3_PKiPKNS_4InfoEiPi+0xa5>
    ecc1:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    ecc5:	48 89 c7             	mov    %rax,%rdi
    ecc8:	e8 03 b9 ff ff       	call   a5d0 <_ZdaPv@plt>
  return newcomm;
    eccd:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    ecd1:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ecd5:	48 89 d6             	mov    %rdx,%rsi
    ecd8:	48 89 c7             	mov    %rax,%rdi
    ecdb:	e8 d2 d0 ff ff       	call   bdb2 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
}
    ece0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ece4:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    eceb:	00 00 
    eced:	74 05                	je     ecf4 <_ZN3MPI9Intracomm14Spawn_multipleEiPPKcPS3_PKiPKNS_4InfoEiPi+0xcc>
    ecef:	e8 ec b7 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    ecf4:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ecf8:	c9                   	leave  
    ecf9:	c3                   	ret    

000000000000ecfa <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t>:
//
//   ========   Cartcomm member functions  ========
//

inline
MPI::Cartcomm::Cartcomm(const MPI_Comm& data) {
    ecfa:	f3 0f 1e fa          	endbr64 
    ecfe:	55                   	push   %rbp
    ecff:	48 89 e5             	mov    %rsp,%rbp
    ed02:	53                   	push   %rbx
    ed03:	48 83 ec 28          	sub    $0x28,%rsp
    ed07:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    ed0b:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    ed0f:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ed16:	00 00 
    ed18:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    ed1c:	31 c0                	xor    %eax,%eax
    ed1e:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ed22:	48 89 c7             	mov    %rax,%rdi
    ed25:	e8 fc cf ff ff       	call   bd26 <_ZN3MPI9IntracommC1Ev>
    ed2a:	48 8d 15 b7 ce 00 00 	lea    0xceb7(%rip),%rdx        # 1bbe8 <_ZTVN3MPI8CartcommE+0x10>
    ed31:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ed35:	48 89 10             	mov    %rdx,(%rax)
  int status = 0;
    ed38:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
  if (MPI::Is_initialized() && (data != MPI_COMM_NULL)) {
    ed3f:	e8 88 db ff ff       	call   c8cc <_ZN3MPI14Is_initializedEv>
    ed44:	84 c0                	test   %al,%al
    ed46:	74 1a                	je     ed62 <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0x68>
    ed48:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ed4c:	48 8b 00             	mov    (%rax),%rax
    ed4f:	48 8d 15 ca f0 00 00 	lea    0xf0ca(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    ed56:	48 39 d0             	cmp    %rdx,%rax
    ed59:	74 07                	je     ed62 <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0x68>
    ed5b:	b8 01 00 00 00       	mov    $0x1,%eax
    ed60:	eb 05                	jmp    ed67 <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0x6d>
    ed62:	b8 00 00 00 00       	mov    $0x0,%eax
    ed67:	84 c0                	test   %al,%al
    ed69:	74 40                	je     edab <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0xb1>
    (void)MPI_Topo_test(data, &status) ;
    ed6b:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ed6f:	48 8b 00             	mov    (%rax),%rax
    ed72:	48 8d 55 e4          	lea    -0x1c(%rbp),%rdx
    ed76:	48 89 d6             	mov    %rdx,%rsi
    ed79:	48 89 c7             	mov    %rax,%rdi
    ed7c:	e8 bf b9 ff ff       	call   a740 <MPI_Topo_test@plt>
    if (status == MPI_CART)
    ed81:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    ed84:	83 f8 01             	cmp    $0x1,%eax
    ed87:	75 11                	jne    ed9a <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0xa0>
      mpi_comm = data;
    ed89:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ed8d:	48 8b 10             	mov    (%rax),%rdx
    ed90:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ed94:	48 89 50 08          	mov    %rdx,0x8(%rax)
      mpi_comm = MPI_COMM_NULL;
  }
  else {
    mpi_comm = data;
  }
}
    ed98:	eb 40                	jmp    edda <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0xe0>
      mpi_comm = MPI_COMM_NULL;
    ed9a:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ed9e:	48 8d 15 7b f0 00 00 	lea    0xf07b(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    eda5:	48 89 50 08          	mov    %rdx,0x8(%rax)
}
    eda9:	eb 2f                	jmp    edda <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0xe0>
    mpi_comm = data;
    edab:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    edaf:	48 8b 10             	mov    (%rax),%rdx
    edb2:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    edb6:	48 89 50 08          	mov    %rdx,0x8(%rax)
}
    edba:	eb 1e                	jmp    edda <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0xe0>
    edbc:	f3 0f 1e fa          	endbr64 
    edc0:	48 89 c3             	mov    %rax,%rbx
    edc3:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    edc7:	48 89 c7             	mov    %rax,%rdi
    edca:	e8 85 cf ff ff       	call   bd54 <_ZN3MPI9IntracommD1Ev>
    edcf:	48 89 d8             	mov    %rbx,%rax
    edd2:	48 89 c7             	mov    %rax,%rdi
    edd5:	e8 d6 bb ff ff       	call   a9b0 <_Unwind_Resume@plt>
    edda:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    edde:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    ede5:	00 00 
    ede7:	74 05                	je     edee <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0xf4>
    ede9:	e8 f2 b6 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    edee:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    edf2:	c9                   	leave  
    edf3:	c3                   	ret    

000000000000edf4 <_ZNK3MPI8Cartcomm7Get_dimEv>:
//  Process Topologies
//

inline int
MPI::Cartcomm::Get_dim() const
{
    edf4:	f3 0f 1e fa          	endbr64 
    edf8:	55                   	push   %rbp
    edf9:	48 89 e5             	mov    %rsp,%rbp
    edfc:	48 83 ec 20          	sub    $0x20,%rsp
    ee00:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    ee04:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ee0b:	00 00 
    ee0d:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    ee11:	31 c0                	xor    %eax,%eax
  int ndims;
  (void)MPI_Cartdim_get(mpi_comm, &ndims);
    ee13:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ee17:	48 8b 40 08          	mov    0x8(%rax),%rax
    ee1b:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    ee1f:	48 89 d6             	mov    %rdx,%rsi
    ee22:	48 89 c7             	mov    %rax,%rdi
    ee25:	e8 16 b6 ff ff       	call   a440 <MPI_Cartdim_get@plt>
  return ndims;
    ee2a:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    ee2d:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    ee31:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    ee38:	00 00 
    ee3a:	74 05                	je     ee41 <_ZNK3MPI8Cartcomm7Get_dimEv+0x4d>
    ee3c:	e8 9f b6 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    ee41:	c9                   	leave  
    ee42:	c3                   	ret    
    ee43:	90                   	nop

000000000000ee44 <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_>:

inline void
MPI::Cartcomm::Get_topo(int maxdims, int dims[], bool periods[],
			       int coords[]) const
{
    ee44:	f3 0f 1e fa          	endbr64 
    ee48:	55                   	push   %rbp
    ee49:	48 89 e5             	mov    %rsp,%rbp
    ee4c:	48 83 ec 40          	sub    $0x40,%rsp
    ee50:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    ee54:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    ee57:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    ee5b:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    ee5f:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
  int *int_periods = new int [maxdims];
    ee63:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    ee66:	48 98                	cltq   
    ee68:	48 ba fe ff ff ff ff 	movabs $0x1ffffffffffffffe,%rdx
    ee6f:	ff ff 1f 
    ee72:	48 39 d0             	cmp    %rdx,%rax
    ee75:	77 19                	ja     ee90 <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_+0x4c>
    ee77:	48 c1 e0 02          	shl    $0x2,%rax
    ee7b:	48 89 c7             	mov    %rax,%rdi
    ee7e:	e8 3d af ff ff       	call   9dc0 <_Znam@plt>
    ee83:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
  int i;
  for (i=0; i<maxdims; i++) {
    ee87:	c7 45 f4 00 00 00 00 	movl   $0x0,-0xc(%rbp)
    ee8e:	eb 33                	jmp    eec3 <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_+0x7f>
  int *int_periods = new int [maxdims];
    ee90:	e8 5b b8 ff ff       	call   a6f0 <__cxa_throw_bad_array_new_length@plt>
    int_periods[i] = (int)periods[i];
    ee95:	8b 45 f4             	mov    -0xc(%rbp),%eax
    ee98:	48 63 d0             	movslq %eax,%rdx
    ee9b:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ee9f:	48 01 d0             	add    %rdx,%rax
    eea2:	0f b6 00             	movzbl (%rax),%eax
    eea5:	8b 55 f4             	mov    -0xc(%rbp),%edx
    eea8:	48 63 d2             	movslq %edx,%rdx
    eeab:	48 8d 0c 95 00 00 00 	lea    0x0(,%rdx,4),%rcx
    eeb2:	00 
    eeb3:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    eeb7:	48 01 ca             	add    %rcx,%rdx
    eeba:	0f b6 c0             	movzbl %al,%eax
    eebd:	89 02                	mov    %eax,(%rdx)
  for (i=0; i<maxdims; i++) {
    eebf:	83 45 f4 01          	addl   $0x1,-0xc(%rbp)
    eec3:	8b 45 f4             	mov    -0xc(%rbp),%eax
    eec6:	3b 45 e4             	cmp    -0x1c(%rbp),%eax
    eec9:	7c ca                	jl     ee95 <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_+0x51>
  }
  (void)MPI_Cart_get(mpi_comm, maxdims, dims, int_periods, coords);
    eecb:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    eecf:	48 8b 40 08          	mov    0x8(%rax),%rax
    eed3:	48 8b 7d c8          	mov    -0x38(%rbp),%rdi
    eed7:	48 8b 4d f8          	mov    -0x8(%rbp),%rcx
    eedb:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
    eedf:	8b 75 e4             	mov    -0x1c(%rbp),%esi
    eee2:	49 89 f8             	mov    %rdi,%r8
    eee5:	48 89 c7             	mov    %rax,%rdi
    eee8:	e8 33 b4 ff ff       	call   a320 <MPI_Cart_get@plt>
  for (i=0; i<maxdims; i++) {
    eeed:	c7 45 f4 00 00 00 00 	movl   $0x0,-0xc(%rbp)
    eef4:	eb 2e                	jmp    ef24 <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_+0xe0>
    periods[i] = OPAL_INT_TO_BOOL(int_periods[i]);
    eef6:	8b 45 f4             	mov    -0xc(%rbp),%eax
    eef9:	48 98                	cltq   
    eefb:	48 8d 14 85 00 00 00 	lea    0x0(,%rax,4),%rdx
    ef02:	00 
    ef03:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ef07:	48 01 d0             	add    %rdx,%rax
    ef0a:	8b 00                	mov    (%rax),%eax
    ef0c:	8b 55 f4             	mov    -0xc(%rbp),%edx
    ef0f:	48 63 ca             	movslq %edx,%rcx
    ef12:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    ef16:	48 01 ca             	add    %rcx,%rdx
    ef19:	85 c0                	test   %eax,%eax
    ef1b:	0f 95 c0             	setne  %al
    ef1e:	88 02                	mov    %al,(%rdx)
  for (i=0; i<maxdims; i++) {
    ef20:	83 45 f4 01          	addl   $0x1,-0xc(%rbp)
    ef24:	8b 45 f4             	mov    -0xc(%rbp),%eax
    ef27:	3b 45 e4             	cmp    -0x1c(%rbp),%eax
    ef2a:	7c ca                	jl     eef6 <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_+0xb2>
  }
  delete [] int_periods;
    ef2c:	48 83 7d f8 00       	cmpq   $0x0,-0x8(%rbp)
    ef31:	74 0c                	je     ef3f <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_+0xfb>
    ef33:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ef37:	48 89 c7             	mov    %rax,%rdi
    ef3a:	e8 91 b6 ff ff       	call   a5d0 <_ZdaPv@plt>
}
    ef3f:	90                   	nop
    ef40:	c9                   	leave  
    ef41:	c3                   	ret    

000000000000ef42 <_ZNK3MPI8Cartcomm13Get_cart_rankEPKi>:

inline int
MPI::Cartcomm::Get_cart_rank(const int coords[]) const
{
    ef42:	f3 0f 1e fa          	endbr64 
    ef46:	55                   	push   %rbp
    ef47:	48 89 e5             	mov    %rsp,%rbp
    ef4a:	48 83 ec 20          	sub    $0x20,%rsp
    ef4e:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    ef52:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    ef56:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ef5d:	00 00 
    ef5f:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    ef63:	31 c0                	xor    %eax,%eax
  int myrank;
  (void)MPI_Cart_rank(mpi_comm, const_cast<int *>(coords), &myrank);
    ef65:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ef69:	48 8b 40 08          	mov    0x8(%rax),%rax
    ef6d:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    ef71:	48 8b 4d e0          	mov    -0x20(%rbp),%rcx
    ef75:	48 89 ce             	mov    %rcx,%rsi
    ef78:	48 89 c7             	mov    %rax,%rdi
    ef7b:	e8 50 ba ff ff       	call   a9d0 <MPI_Cart_rank@plt>
  return myrank;
    ef80:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    ef83:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    ef87:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    ef8e:	00 00 
    ef90:	74 05                	je     ef97 <_ZNK3MPI8Cartcomm13Get_cart_rankEPKi+0x55>
    ef92:	e8 49 b5 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    ef97:	c9                   	leave  
    ef98:	c3                   	ret    
    ef99:	90                   	nop

000000000000ef9a <_ZNK3MPI8Cartcomm10Get_coordsEiiPi>:

inline void
MPI::Cartcomm::Get_coords(int rank, int maxdims, int coords[]) const
{
    ef9a:	f3 0f 1e fa          	endbr64 
    ef9e:	55                   	push   %rbp
    ef9f:	48 89 e5             	mov    %rsp,%rbp
    efa2:	48 83 ec 20          	sub    $0x20,%rsp
    efa6:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    efaa:	89 75 f4             	mov    %esi,-0xc(%rbp)
    efad:	89 55 f0             	mov    %edx,-0x10(%rbp)
    efb0:	48 89 4d e8          	mov    %rcx,-0x18(%rbp)
  (void)MPI_Cart_coords(mpi_comm, rank, maxdims, coords);
    efb4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    efb8:	48 8b 40 08          	mov    0x8(%rax),%rax
    efbc:	48 8b 4d e8          	mov    -0x18(%rbp),%rcx
    efc0:	8b 55 f0             	mov    -0x10(%rbp),%edx
    efc3:	8b 75 f4             	mov    -0xc(%rbp),%esi
    efc6:	48 89 c7             	mov    %rax,%rdi
    efc9:	e8 82 b4 ff ff       	call   a450 <MPI_Cart_coords@plt>
}
    efce:	90                   	nop
    efcf:	c9                   	leave  
    efd0:	c3                   	ret    
    efd1:	90                   	nop

000000000000efd2 <_ZNK3MPI8Cartcomm5ShiftEiiRiS1_>:

inline void
MPI::Cartcomm::Shift(int direction, int disp,
			    int &rank_source, int &rank_dest) const
{
    efd2:	f3 0f 1e fa          	endbr64 
    efd6:	55                   	push   %rbp
    efd7:	48 89 e5             	mov    %rsp,%rbp
    efda:	48 83 ec 20          	sub    $0x20,%rsp
    efde:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    efe2:	89 75 f4             	mov    %esi,-0xc(%rbp)
    efe5:	89 55 f0             	mov    %edx,-0x10(%rbp)
    efe8:	48 89 4d e8          	mov    %rcx,-0x18(%rbp)
    efec:	4c 89 45 e0          	mov    %r8,-0x20(%rbp)
  (void)MPI_Cart_shift(mpi_comm, direction, disp, &rank_source, &rank_dest);
    eff0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    eff4:	48 8b 40 08          	mov    0x8(%rax),%rax
    eff8:	48 8b 7d e0          	mov    -0x20(%rbp),%rdi
    effc:	48 8b 4d e8          	mov    -0x18(%rbp),%rcx
    f000:	8b 55 f0             	mov    -0x10(%rbp),%edx
    f003:	8b 75 f4             	mov    -0xc(%rbp),%esi
    f006:	49 89 f8             	mov    %rdi,%r8
    f009:	48 89 c7             	mov    %rax,%rdi
    f00c:	e8 cf b0 ff ff       	call   a0e0 <MPI_Cart_shift@plt>
}
    f011:	90                   	nop
    f012:	c9                   	leave  
    f013:	c3                   	ret    

000000000000f014 <_ZNK3MPI8Cartcomm3SubEPKb>:

inline MPI::Cartcomm
MPI::Cartcomm::Sub(const bool remain_dims[]) const
{
    f014:	f3 0f 1e fa          	endbr64 
    f018:	55                   	push   %rbp
    f019:	48 89 e5             	mov    %rsp,%rbp
    f01c:	48 83 ec 40          	sub    $0x40,%rsp
    f020:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    f024:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    f028:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    f02c:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f033:	00 00 
    f035:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f039:	31 c0                	xor    %eax,%eax
  int ndims;
  MPI_Cartdim_get(mpi_comm, &ndims);
    f03b:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    f03f:	48 8b 40 08          	mov    0x8(%rax),%rax
    f043:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    f047:	48 89 d6             	mov    %rdx,%rsi
    f04a:	48 89 c7             	mov    %rax,%rdi
    f04d:	e8 ee b3 ff ff       	call   a440 <MPI_Cartdim_get@plt>
  int* int_remain_dims = new int[ndims];
    f052:	8b 45 e0             	mov    -0x20(%rbp),%eax
    f055:	48 98                	cltq   
    f057:	48 ba fe ff ff ff ff 	movabs $0x1ffffffffffffffe,%rdx
    f05e:	ff ff 1f 
    f061:	48 39 d0             	cmp    %rdx,%rax
    f064:	77 19                	ja     f07f <_ZNK3MPI8Cartcomm3SubEPKb+0x6b>
    f066:	48 c1 e0 02          	shl    $0x2,%rax
    f06a:	48 89 c7             	mov    %rax,%rdi
    f06d:	e8 4e ad ff ff       	call   9dc0 <_Znam@plt>
    f072:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
  for (int i=0; i<ndims; i++) {
    f076:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
    f07d:	eb 33                	jmp    f0b2 <_ZNK3MPI8Cartcomm3SubEPKb+0x9e>
  int* int_remain_dims = new int[ndims];
    f07f:	e8 6c b6 ff ff       	call   a6f0 <__cxa_throw_bad_array_new_length@plt>
    int_remain_dims[i] = (int)remain_dims[i];
    f084:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    f087:	48 63 d0             	movslq %eax,%rdx
    f08a:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    f08e:	48 01 d0             	add    %rdx,%rax
    f091:	0f b6 00             	movzbl (%rax),%eax
    f094:	8b 55 e4             	mov    -0x1c(%rbp),%edx
    f097:	48 63 d2             	movslq %edx,%rdx
    f09a:	48 8d 0c 95 00 00 00 	lea    0x0(,%rdx,4),%rcx
    f0a1:	00 
    f0a2:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f0a6:	48 01 ca             	add    %rcx,%rdx
    f0a9:	0f b6 c0             	movzbl %al,%eax
    f0ac:	89 02                	mov    %eax,(%rdx)
  for (int i=0; i<ndims; i++) {
    f0ae:	83 45 e4 01          	addl   $0x1,-0x1c(%rbp)
    f0b2:	8b 45 e0             	mov    -0x20(%rbp),%eax
    f0b5:	39 45 e4             	cmp    %eax,-0x1c(%rbp)
    f0b8:	7c ca                	jl     f084 <_ZNK3MPI8Cartcomm3SubEPKb+0x70>
  }
  MPI_Comm newcomm;
  (void)MPI_Cart_sub(mpi_comm, int_remain_dims, &newcomm);
    f0ba:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    f0be:	48 8b 40 08          	mov    0x8(%rax),%rax
    f0c2:	48 8d 55 e8          	lea    -0x18(%rbp),%rdx
    f0c6:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    f0ca:	48 89 ce             	mov    %rcx,%rsi
    f0cd:	48 89 c7             	mov    %rax,%rdi
    f0d0:	e8 6b b0 ff ff       	call   a140 <MPI_Cart_sub@plt>
  delete [] int_remain_dims;
    f0d5:	48 83 7d f0 00       	cmpq   $0x0,-0x10(%rbp)
    f0da:	74 0c                	je     f0e8 <_ZNK3MPI8Cartcomm3SubEPKb+0xd4>
    f0dc:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    f0e0:	48 89 c7             	mov    %rax,%rdi
    f0e3:	e8 e8 b4 ff ff       	call   a5d0 <_ZdaPv@plt>
  return newcomm;
    f0e8:	48 8d 55 e8          	lea    -0x18(%rbp),%rdx
    f0ec:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f0f0:	48 89 d6             	mov    %rdx,%rsi
    f0f3:	48 89 c7             	mov    %rax,%rdi
    f0f6:	e8 ff fb ff ff       	call   ecfa <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t>
}
    f0fb:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f0ff:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f106:	00 00 
    f108:	74 05                	je     f10f <_ZNK3MPI8Cartcomm3SubEPKb+0xfb>
    f10a:	e8 d1 b3 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f10f:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f113:	c9                   	leave  
    f114:	c3                   	ret    
    f115:	90                   	nop

000000000000f116 <_ZNK3MPI8Cartcomm3MapEiPKiPKb>:

inline int
MPI::Cartcomm::Map(int ndims, const int dims[], const bool periods[]) const
{
    f116:	f3 0f 1e fa          	endbr64 
    f11a:	55                   	push   %rbp
    f11b:	48 89 e5             	mov    %rsp,%rbp
    f11e:	48 83 ec 40          	sub    $0x40,%rsp
    f122:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    f126:	89 75 d4             	mov    %esi,-0x2c(%rbp)
    f129:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    f12d:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    f131:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f138:	00 00 
    f13a:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f13e:	31 c0                	xor    %eax,%eax
  int *int_periods = new int [ndims];
    f140:	8b 45 d4             	mov    -0x2c(%rbp),%eax
    f143:	48 98                	cltq   
    f145:	48 ba fe ff ff ff ff 	movabs $0x1ffffffffffffffe,%rdx
    f14c:	ff ff 1f 
    f14f:	48 39 d0             	cmp    %rdx,%rax
    f152:	77 19                	ja     f16d <_ZNK3MPI8Cartcomm3MapEiPKiPKb+0x57>
    f154:	48 c1 e0 02          	shl    $0x2,%rax
    f158:	48 89 c7             	mov    %rax,%rdi
    f15b:	e8 60 ac ff ff       	call   9dc0 <_Znam@plt>
    f160:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
  for (int i=0; i<ndims; i++) {
    f164:	c7 45 ec 00 00 00 00 	movl   $0x0,-0x14(%rbp)
    f16b:	eb 33                	jmp    f1a0 <_ZNK3MPI8Cartcomm3MapEiPKiPKb+0x8a>
  int *int_periods = new int [ndims];
    f16d:	e8 7e b5 ff ff       	call   a6f0 <__cxa_throw_bad_array_new_length@plt>
    int_periods[i] = (int) periods[i];
    f172:	8b 45 ec             	mov    -0x14(%rbp),%eax
    f175:	48 63 d0             	movslq %eax,%rdx
    f178:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    f17c:	48 01 d0             	add    %rdx,%rax
    f17f:	0f b6 00             	movzbl (%rax),%eax
    f182:	8b 55 ec             	mov    -0x14(%rbp),%edx
    f185:	48 63 d2             	movslq %edx,%rdx
    f188:	48 8d 0c 95 00 00 00 	lea    0x0(,%rdx,4),%rcx
    f18f:	00 
    f190:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f194:	48 01 ca             	add    %rcx,%rdx
    f197:	0f b6 c0             	movzbl %al,%eax
    f19a:	89 02                	mov    %eax,(%rdx)
  for (int i=0; i<ndims; i++) {
    f19c:	83 45 ec 01          	addl   $0x1,-0x14(%rbp)
    f1a0:	8b 45 ec             	mov    -0x14(%rbp),%eax
    f1a3:	3b 45 d4             	cmp    -0x2c(%rbp),%eax
    f1a6:	7c ca                	jl     f172 <_ZNK3MPI8Cartcomm3MapEiPKiPKb+0x5c>
  }
  int newrank;
  (void)MPI_Cart_map(mpi_comm, ndims, const_cast<int *>(dims), int_periods, &newrank);
    f1a8:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f1ac:	48 8b 40 08          	mov    0x8(%rax),%rax
    f1b0:	48 8d 7d e8          	lea    -0x18(%rbp),%rdi
    f1b4:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    f1b8:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    f1bc:	8b 75 d4             	mov    -0x2c(%rbp),%esi
    f1bf:	49 89 f8             	mov    %rdi,%r8
    f1c2:	48 89 c7             	mov    %rax,%rdi
    f1c5:	e8 c6 b1 ff ff       	call   a390 <MPI_Cart_map@plt>
  delete [] int_periods;
    f1ca:	48 83 7d f0 00       	cmpq   $0x0,-0x10(%rbp)
    f1cf:	74 0c                	je     f1dd <_ZNK3MPI8Cartcomm3MapEiPKiPKb+0xc7>
    f1d1:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    f1d5:	48 89 c7             	mov    %rax,%rdi
    f1d8:	e8 f3 b3 ff ff       	call   a5d0 <_ZdaPv@plt>
  return newrank;
    f1dd:	8b 45 e8             	mov    -0x18(%rbp),%eax
}
    f1e0:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    f1e4:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f1eb:	00 00 
    f1ed:	74 05                	je     f1f4 <_ZNK3MPI8Cartcomm3MapEiPKiPKb+0xde>
    f1ef:	e8 ec b2 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f1f4:	c9                   	leave  
    f1f5:	c3                   	ret    

000000000000f1f6 <_ZNK3MPI8Cartcomm5CloneEv>:


inline MPI::Cartcomm&
MPI::Cartcomm::Clone() const
{
    f1f6:	f3 0f 1e fa          	endbr64 
    f1fa:	55                   	push   %rbp
    f1fb:	48 89 e5             	mov    %rsp,%rbp
    f1fe:	41 54                	push   %r12
    f200:	53                   	push   %rbx
    f201:	48 83 ec 30          	sub    $0x30,%rsp
    f205:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
    f209:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f210:	00 00 
    f212:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    f216:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void)MPI_Comm_dup(mpi_comm, &newcomm);
    f218:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    f21c:	48 8b 40 08          	mov    0x8(%rax),%rax
    f220:	48 8d 55 d8          	lea    -0x28(%rbp),%rdx
    f224:	48 89 d6             	mov    %rdx,%rsi
    f227:	48 89 c7             	mov    %rax,%rdi
    f22a:	e8 21 b4 ff ff       	call   a650 <MPI_Comm_dup@plt>
  MPI::Cartcomm* dup = new MPI::Cartcomm(newcomm);
    f22f:	bf 10 00 00 00       	mov    $0x10,%edi
    f234:	e8 a7 b1 ff ff       	call   a3e0 <_Znwm@plt>
    f239:	48 89 c3             	mov    %rax,%rbx
    f23c:	48 8d 45 d8          	lea    -0x28(%rbp),%rax
    f240:	48 89 c6             	mov    %rax,%rsi
    f243:	48 89 df             	mov    %rbx,%rdi
    f246:	e8 af fa ff ff       	call   ecfa <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t>
    f24b:	48 89 5d e0          	mov    %rbx,-0x20(%rbp)
  return *dup;
    f24f:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
}
    f253:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    f257:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f25e:	00 00 
    f260:	74 26                	je     f288 <_ZNK3MPI8Cartcomm5CloneEv+0x92>
    f262:	eb 1f                	jmp    f283 <_ZNK3MPI8Cartcomm5CloneEv+0x8d>
    f264:	f3 0f 1e fa          	endbr64 
  MPI::Cartcomm* dup = new MPI::Cartcomm(newcomm);
    f268:	49 89 c4             	mov    %rax,%r12
    f26b:	be 10 00 00 00       	mov    $0x10,%esi
    f270:	48 89 df             	mov    %rbx,%rdi
    f273:	e8 88 b1 ff ff       	call   a400 <_ZdlPvm@plt>
    f278:	4c 89 e0             	mov    %r12,%rax
    f27b:	48 89 c7             	mov    %rax,%rdi
    f27e:	e8 2d b7 ff ff       	call   a9b0 <_Unwind_Resume@plt>
}
    f283:	e8 58 b2 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f288:	48 83 c4 30          	add    $0x30,%rsp
    f28c:	5b                   	pop    %rbx
    f28d:	41 5c                	pop    %r12
    f28f:	5d                   	pop    %rbp
    f290:	c3                   	ret    
    f291:	90                   	nop

000000000000f292 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t>:
//
//   ========   Graphcomm member functions  ========
//

inline
MPI::Graphcomm::Graphcomm(const MPI_Comm& data) {
    f292:	f3 0f 1e fa          	endbr64 
    f296:	55                   	push   %rbp
    f297:	48 89 e5             	mov    %rsp,%rbp
    f29a:	53                   	push   %rbx
    f29b:	48 83 ec 28          	sub    $0x28,%rsp
    f29f:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    f2a3:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    f2a7:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f2ae:	00 00 
    f2b0:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    f2b4:	31 c0                	xor    %eax,%eax
    f2b6:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f2ba:	48 89 c7             	mov    %rax,%rdi
    f2bd:	e8 64 ca ff ff       	call   bd26 <_ZN3MPI9IntracommC1Ev>
    f2c2:	48 8d 15 bf c6 00 00 	lea    0xc6bf(%rip),%rdx        # 1b988 <_ZTVN3MPI9GraphcommE+0x10>
    f2c9:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f2cd:	48 89 10             	mov    %rdx,(%rax)
  int status = 0;
    f2d0:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
  if (MPI::Is_initialized() && (data != MPI_COMM_NULL)) {
    f2d7:	e8 f0 d5 ff ff       	call   c8cc <_ZN3MPI14Is_initializedEv>
    f2dc:	84 c0                	test   %al,%al
    f2de:	74 1a                	je     f2fa <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0x68>
    f2e0:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    f2e4:	48 8b 00             	mov    (%rax),%rax
    f2e7:	48 8d 15 32 eb 00 00 	lea    0xeb32(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    f2ee:	48 39 d0             	cmp    %rdx,%rax
    f2f1:	74 07                	je     f2fa <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0x68>
    f2f3:	b8 01 00 00 00       	mov    $0x1,%eax
    f2f8:	eb 05                	jmp    f2ff <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0x6d>
    f2fa:	b8 00 00 00 00       	mov    $0x0,%eax
    f2ff:	84 c0                	test   %al,%al
    f301:	74 40                	je     f343 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0xb1>
    (void)MPI_Topo_test(data, &status) ;
    f303:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    f307:	48 8b 00             	mov    (%rax),%rax
    f30a:	48 8d 55 e4          	lea    -0x1c(%rbp),%rdx
    f30e:	48 89 d6             	mov    %rdx,%rsi
    f311:	48 89 c7             	mov    %rax,%rdi
    f314:	e8 27 b4 ff ff       	call   a740 <MPI_Topo_test@plt>
    if (status == MPI_GRAPH)
    f319:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    f31c:	83 f8 02             	cmp    $0x2,%eax
    f31f:	75 11                	jne    f332 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0xa0>
      mpi_comm = data;
    f321:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    f325:	48 8b 10             	mov    (%rax),%rdx
    f328:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f32c:	48 89 50 08          	mov    %rdx,0x8(%rax)
      mpi_comm = MPI_COMM_NULL;
  }
  else {
    mpi_comm = data;
  }
}
    f330:	eb 40                	jmp    f372 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0xe0>
      mpi_comm = MPI_COMM_NULL;
    f332:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f336:	48 8d 15 e3 ea 00 00 	lea    0xeae3(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    f33d:	48 89 50 08          	mov    %rdx,0x8(%rax)
}
    f341:	eb 2f                	jmp    f372 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0xe0>
    mpi_comm = data;
    f343:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    f347:	48 8b 10             	mov    (%rax),%rdx
    f34a:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f34e:	48 89 50 08          	mov    %rdx,0x8(%rax)
}
    f352:	eb 1e                	jmp    f372 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0xe0>
    f354:	f3 0f 1e fa          	endbr64 
    f358:	48 89 c3             	mov    %rax,%rbx
    f35b:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f35f:	48 89 c7             	mov    %rax,%rdi
    f362:	e8 ed c9 ff ff       	call   bd54 <_ZN3MPI9IntracommD1Ev>
    f367:	48 89 d8             	mov    %rbx,%rax
    f36a:	48 89 c7             	mov    %rax,%rdi
    f36d:	e8 3e b6 ff ff       	call   a9b0 <_Unwind_Resume@plt>
    f372:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f376:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f37d:	00 00 
    f37f:	74 05                	je     f386 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0xf4>
    f381:	e8 5a b1 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f386:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    f38a:	c9                   	leave  
    f38b:	c3                   	ret    

000000000000f38c <_ZNK3MPI9Graphcomm8Get_dimsEPiS1_>:
//  Process Topologies
//

inline void
MPI::Graphcomm::Get_dims(int nnodes[], int nedges[]) const
{
    f38c:	f3 0f 1e fa          	endbr64 
    f390:	55                   	push   %rbp
    f391:	48 89 e5             	mov    %rsp,%rbp
    f394:	48 83 ec 20          	sub    $0x20,%rsp
    f398:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    f39c:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    f3a0:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
  (void)MPI_Graphdims_get(mpi_comm, nnodes, nedges);
    f3a4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f3a8:	48 8b 40 08          	mov    0x8(%rax),%rax
    f3ac:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    f3b0:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    f3b4:	48 89 ce             	mov    %rcx,%rsi
    f3b7:	48 89 c7             	mov    %rax,%rdi
    f3ba:	e8 01 ad ff ff       	call   a0c0 <MPI_Graphdims_get@plt>
}
    f3bf:	90                   	nop
    f3c0:	c9                   	leave  
    f3c1:	c3                   	ret    

000000000000f3c2 <_ZNK3MPI9Graphcomm8Get_topoEiiPiS1_>:

inline void
MPI::Graphcomm::Get_topo(int maxindex, int maxedges, int index[],
	 int edges[]) const
{
    f3c2:	f3 0f 1e fa          	endbr64 
    f3c6:	55                   	push   %rbp
    f3c7:	48 89 e5             	mov    %rsp,%rbp
    f3ca:	48 83 ec 20          	sub    $0x20,%rsp
    f3ce:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    f3d2:	89 75 f4             	mov    %esi,-0xc(%rbp)
    f3d5:	89 55 f0             	mov    %edx,-0x10(%rbp)
    f3d8:	48 89 4d e8          	mov    %rcx,-0x18(%rbp)
    f3dc:	4c 89 45 e0          	mov    %r8,-0x20(%rbp)
  (void)MPI_Graph_get(mpi_comm, maxindex, maxedges, index, edges);
    f3e0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f3e4:	48 8b 40 08          	mov    0x8(%rax),%rax
    f3e8:	48 8b 7d e0          	mov    -0x20(%rbp),%rdi
    f3ec:	48 8b 4d e8          	mov    -0x18(%rbp),%rcx
    f3f0:	8b 55 f0             	mov    -0x10(%rbp),%edx
    f3f3:	8b 75 f4             	mov    -0xc(%rbp),%esi
    f3f6:	49 89 f8             	mov    %rdi,%r8
    f3f9:	48 89 c7             	mov    %rax,%rdi
    f3fc:	e8 ef ae ff ff       	call   a2f0 <MPI_Graph_get@plt>
}
    f401:	90                   	nop
    f402:	c9                   	leave  
    f403:	c3                   	ret    

000000000000f404 <_ZNK3MPI9Graphcomm19Get_neighbors_countEi>:

inline int
MPI::Graphcomm::Get_neighbors_count(int rank) const
{
    f404:	f3 0f 1e fa          	endbr64 
    f408:	55                   	push   %rbp
    f409:	48 89 e5             	mov    %rsp,%rbp
    f40c:	48 83 ec 20          	sub    $0x20,%rsp
    f410:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f414:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    f417:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f41e:	00 00 
    f420:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f424:	31 c0                	xor    %eax,%eax
  int nneighbors;
  (void)MPI_Graph_neighbors_count(mpi_comm, rank, &nneighbors);
    f426:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f42a:	48 8b 40 08          	mov    0x8(%rax),%rax
    f42e:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    f432:	8b 4d e4             	mov    -0x1c(%rbp),%ecx
    f435:	89 ce                	mov    %ecx,%esi
    f437:	48 89 c7             	mov    %rax,%rdi
    f43a:	e8 71 b6 ff ff       	call   aab0 <MPI_Graph_neighbors_count@plt>
  return nneighbors;
    f43f:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    f442:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    f446:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f44d:	00 00 
    f44f:	74 05                	je     f456 <_ZNK3MPI9Graphcomm19Get_neighbors_countEi+0x52>
    f451:	e8 8a b0 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f456:	c9                   	leave  
    f457:	c3                   	ret    

000000000000f458 <_ZNK3MPI9Graphcomm13Get_neighborsEiiPi>:

inline void
MPI::Graphcomm::Get_neighbors(int rank, int maxneighbors,
	      int neighbors[]) const
{
    f458:	f3 0f 1e fa          	endbr64 
    f45c:	55                   	push   %rbp
    f45d:	48 89 e5             	mov    %rsp,%rbp
    f460:	48 83 ec 20          	sub    $0x20,%rsp
    f464:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    f468:	89 75 f4             	mov    %esi,-0xc(%rbp)
    f46b:	89 55 f0             	mov    %edx,-0x10(%rbp)
    f46e:	48 89 4d e8          	mov    %rcx,-0x18(%rbp)
  (void)MPI_Graph_neighbors(mpi_comm, rank, maxneighbors, neighbors);
    f472:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f476:	48 8b 40 08          	mov    0x8(%rax),%rax
    f47a:	48 8b 4d e8          	mov    -0x18(%rbp),%rcx
    f47e:	8b 55 f0             	mov    -0x10(%rbp),%edx
    f481:	8b 75 f4             	mov    -0xc(%rbp),%esi
    f484:	48 89 c7             	mov    %rax,%rdi
    f487:	e8 c4 aa ff ff       	call   9f50 <MPI_Graph_neighbors@plt>
}
    f48c:	90                   	nop
    f48d:	c9                   	leave  
    f48e:	c3                   	ret    
    f48f:	90                   	nop

000000000000f490 <_ZNK3MPI9Graphcomm3MapEiPKiS2_>:

inline int
MPI::Graphcomm::Map(int nnodes, const int index[],
    const int edges[]) const
{
    f490:	f3 0f 1e fa          	endbr64 
    f494:	55                   	push   %rbp
    f495:	48 89 e5             	mov    %rsp,%rbp
    f498:	48 83 ec 30          	sub    $0x30,%rsp
    f49c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f4a0:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    f4a3:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    f4a7:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    f4ab:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f4b2:	00 00 
    f4b4:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f4b8:	31 c0                	xor    %eax,%eax
  int newrank;
  (void)MPI_Graph_map(mpi_comm, nnodes, const_cast<int *>(index), const_cast<int *>(edges), &newrank);
    f4ba:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f4be:	48 8b 40 08          	mov    0x8(%rax),%rax
    f4c2:	48 8d 7d f4          	lea    -0xc(%rbp),%rdi
    f4c6:	48 8b 4d d0          	mov    -0x30(%rbp),%rcx
    f4ca:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
    f4ce:	8b 75 e4             	mov    -0x1c(%rbp),%esi
    f4d1:	49 89 f8             	mov    %rdi,%r8
    f4d4:	48 89 c7             	mov    %rax,%rdi
    f4d7:	e8 84 ae ff ff       	call   a360 <MPI_Graph_map@plt>
  return newrank;
    f4dc:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    f4df:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    f4e3:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f4ea:	00 00 
    f4ec:	74 05                	je     f4f3 <_ZNK3MPI9Graphcomm3MapEiPKiS2_+0x63>
    f4ee:	e8 ed af ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f4f3:	c9                   	leave  
    f4f4:	c3                   	ret    
    f4f5:	90                   	nop

000000000000f4f6 <_ZNK3MPI9Graphcomm5CloneEv>:

inline MPI::Graphcomm&
MPI::Graphcomm::Clone() const
{
    f4f6:	f3 0f 1e fa          	endbr64 
    f4fa:	55                   	push   %rbp
    f4fb:	48 89 e5             	mov    %rsp,%rbp
    f4fe:	41 54                	push   %r12
    f500:	53                   	push   %rbx
    f501:	48 83 ec 30          	sub    $0x30,%rsp
    f505:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
    f509:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f510:	00 00 
    f512:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    f516:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void)MPI_Comm_dup(mpi_comm, &newcomm);
    f518:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    f51c:	48 8b 40 08          	mov    0x8(%rax),%rax
    f520:	48 8d 55 d8          	lea    -0x28(%rbp),%rdx
    f524:	48 89 d6             	mov    %rdx,%rsi
    f527:	48 89 c7             	mov    %rax,%rdi
    f52a:	e8 21 b1 ff ff       	call   a650 <MPI_Comm_dup@plt>
  MPI::Graphcomm* dup = new MPI::Graphcomm(newcomm);
    f52f:	bf 10 00 00 00       	mov    $0x10,%edi
    f534:	e8 a7 ae ff ff       	call   a3e0 <_Znwm@plt>
    f539:	48 89 c3             	mov    %rax,%rbx
    f53c:	48 8d 45 d8          	lea    -0x28(%rbp),%rax
    f540:	48 89 c6             	mov    %rax,%rsi
    f543:	48 89 df             	mov    %rbx,%rdi
    f546:	e8 47 fd ff ff       	call   f292 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t>
    f54b:	48 89 5d e0          	mov    %rbx,-0x20(%rbp)
  return *dup;
    f54f:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
}
    f553:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    f557:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f55e:	00 00 
    f560:	74 26                	je     f588 <_ZNK3MPI9Graphcomm5CloneEv+0x92>
    f562:	eb 1f                	jmp    f583 <_ZNK3MPI9Graphcomm5CloneEv+0x8d>
    f564:	f3 0f 1e fa          	endbr64 
  MPI::Graphcomm* dup = new MPI::Graphcomm(newcomm);
    f568:	49 89 c4             	mov    %rax,%r12
    f56b:	be 10 00 00 00       	mov    $0x10,%esi
    f570:	48 89 df             	mov    %rbx,%rdi
    f573:	e8 88 ae ff ff       	call   a400 <_ZdlPvm@plt>
    f578:	4c 89 e0             	mov    %r12,%rax
    f57b:	48 89 c7             	mov    %rax,%rdi
    f57e:	e8 2d b4 ff ff       	call   a9b0 <_Unwind_Resume@plt>
}
    f583:	e8 58 af ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f588:	48 83 c4 30          	add    $0x30,%rsp
    f58c:	5b                   	pop    %rbx
    f58d:	41 5c                	pop    %r12
    f58f:	5d                   	pop    %rbp
    f590:	c3                   	ret    
    f591:	90                   	nop

000000000000f592 <_ZNK3MPI9Intercomm5CloneEv>:
  return newcomm;
}

inline MPI::Intercomm&
MPI::Intercomm::Clone() const
{
    f592:	f3 0f 1e fa          	endbr64 
    f596:	55                   	push   %rbp
    f597:	48 89 e5             	mov    %rsp,%rbp
    f59a:	53                   	push   %rbx
    f59b:	48 83 ec 38          	sub    $0x38,%rsp
    f59f:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
    f5a3:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f5aa:	00 00 
    f5ac:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    f5b0:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void)MPI_Comm_dup(mpi_comm, &newcomm);
    f5b2:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    f5b6:	48 8b 40 08          	mov    0x8(%rax),%rax
    f5ba:	48 8d 55 d8          	lea    -0x28(%rbp),%rdx
    f5be:	48 89 d6             	mov    %rdx,%rsi
    f5c1:	48 89 c7             	mov    %rax,%rdi
    f5c4:	e8 87 b0 ff ff       	call   a650 <MPI_Comm_dup@plt>
  MPI::Intercomm* dup = new MPI::Intercomm(newcomm);
    f5c9:	bf 10 00 00 00       	mov    $0x10,%edi
    f5ce:	e8 0d ae ff ff       	call   a3e0 <_Znwm@plt>
    f5d3:	48 89 c3             	mov    %rax,%rbx
    f5d6:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f5da:	48 89 c6             	mov    %rax,%rsi
    f5dd:	48 89 df             	mov    %rbx,%rdi
    f5e0:	e8 cd c7 ff ff       	call   bdb2 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
    f5e5:	48 89 5d e0          	mov    %rbx,-0x20(%rbp)
  return *dup;
    f5e9:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
}
    f5ed:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    f5f1:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f5f8:	00 00 
    f5fa:	74 05                	je     f601 <_ZNK3MPI9Intercomm5CloneEv+0x6f>
    f5fc:	e8 df ae ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f601:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    f605:	c9                   	leave  
    f606:	c3                   	ret    
    f607:	90                   	nop

000000000000f608 <_ZNK3MPI9Intercomm15Get_remote_sizeEv>:

inline int
MPI::Intercomm::Get_remote_size() const
{
    f608:	f3 0f 1e fa          	endbr64 
    f60c:	55                   	push   %rbp
    f60d:	48 89 e5             	mov    %rsp,%rbp
    f610:	48 83 ec 20          	sub    $0x20,%rsp
    f614:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f618:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f61f:	00 00 
    f621:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f625:	31 c0                	xor    %eax,%eax
  int size;
  (void)MPI_Comm_remote_size(mpi_comm, &size);
    f627:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f62b:	48 8b 40 08          	mov    0x8(%rax),%rax
    f62f:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    f633:	48 89 d6             	mov    %rdx,%rsi
    f636:	48 89 c7             	mov    %rax,%rdi
    f639:	e8 62 b1 ff ff       	call   a7a0 <MPI_Comm_remote_size@plt>
  return size;
    f63e:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    f641:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    f645:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f64c:	00 00 
    f64e:	74 05                	je     f655 <_ZNK3MPI9Intercomm15Get_remote_sizeEv+0x4d>
    f650:	e8 8b ae ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f655:	c9                   	leave  
    f656:	c3                   	ret    
    f657:	90                   	nop

000000000000f658 <_ZNK3MPI9Intercomm16Get_remote_groupEv>:

inline MPI::Group
MPI::Intercomm::Get_remote_group() const
{
    f658:	f3 0f 1e fa          	endbr64 
    f65c:	55                   	push   %rbp
    f65d:	48 89 e5             	mov    %rsp,%rbp
    f660:	48 83 ec 20          	sub    $0x20,%rsp
    f664:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f668:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f66c:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f673:	00 00 
    f675:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f679:	31 c0                	xor    %eax,%eax
  MPI_Group group;
  (void)MPI_Comm_remote_group(mpi_comm, &group);
    f67b:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f67f:	48 8b 40 08          	mov    0x8(%rax),%rax
    f683:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    f687:	48 89 d6             	mov    %rdx,%rsi
    f68a:	48 89 c7             	mov    %rax,%rdi
    f68d:	e8 8e af ff ff       	call   a620 <MPI_Comm_remote_group@plt>
  return group;
    f692:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f696:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f69a:	48 89 d6             	mov    %rdx,%rsi
    f69d:	48 89 c7             	mov    %rax,%rdi
    f6a0:	e8 e7 c3 ff ff       	call   ba8c <_ZN3MPI5GroupC1EP12ompi_group_t>
}
    f6a5:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f6a9:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f6b0:	00 00 
    f6b2:	74 05                	je     f6b9 <_ZNK3MPI9Intercomm16Get_remote_groupEv+0x61>
    f6b4:	e8 27 ae ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f6b9:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f6bd:	c9                   	leave  
    f6be:	c3                   	ret    
    f6bf:	90                   	nop

000000000000f6c0 <_ZNK3MPI9Intercomm5MergeEb>:

inline MPI::Intracomm
MPI::Intercomm::Merge(bool high) const
{
    f6c0:	f3 0f 1e fa          	endbr64 
    f6c4:	55                   	push   %rbp
    f6c5:	48 89 e5             	mov    %rsp,%rbp
    f6c8:	48 83 ec 30          	sub    $0x30,%rsp
    f6cc:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f6d0:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f6d4:	89 d0                	mov    %edx,%eax
    f6d6:	88 45 dc             	mov    %al,-0x24(%rbp)
    f6d9:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f6e0:	00 00 
    f6e2:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f6e6:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void)MPI_Intercomm_merge(mpi_comm, (int)high, &newcomm);
    f6e8:	0f b6 4d dc          	movzbl -0x24(%rbp),%ecx
    f6ec:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f6f0:	48 8b 40 08          	mov    0x8(%rax),%rax
    f6f4:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    f6f8:	89 ce                	mov    %ecx,%esi
    f6fa:	48 89 c7             	mov    %rax,%rdi
    f6fd:	e8 4e b3 ff ff       	call   aa50 <MPI_Intercomm_merge@plt>
  return newcomm;
    f702:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f706:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f70a:	48 89 d6             	mov    %rdx,%rsi
    f70d:	48 89 c7             	mov    %rax,%rdi
    f710:	e8 cd eb ff ff       	call   e2e2 <_ZN3MPI9IntracommC1EP19ompi_communicator_t>
}
    f715:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f719:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f720:	00 00 
    f722:	74 05                	je     f729 <_ZNK3MPI9Intercomm5MergeEb+0x69>
    f724:	e8 b7 ad ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f729:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f72d:	c9                   	leave  
    f72e:	c3                   	ret    
    f72f:	90                   	nop

000000000000f730 <_ZNK3MPI9Intercomm6CreateERKNS_5GroupE>:
// Extended Collective Operations
//

inline MPI::Intercomm
MPI::Intercomm::Create(const Group& group) const
{
    f730:	f3 0f 1e fa          	endbr64 
    f734:	55                   	push   %rbp
    f735:	48 89 e5             	mov    %rsp,%rbp
    f738:	48 83 ec 30          	sub    $0x30,%rsp
    f73c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f740:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f744:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    f748:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f74f:	00 00 
    f751:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f755:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void) MPI_Comm_create(mpi_comm, (MPI_Group) group, &newcomm);
    f757:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f75b:	48 89 c7             	mov    %rax,%rdi
    f75e:	e8 a5 c3 ff ff       	call   bb08 <_ZNK3MPI5GroupcvP12ompi_group_tEv>
    f763:	48 89 c1             	mov    %rax,%rcx
    f766:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f76a:	48 8b 40 08          	mov    0x8(%rax),%rax
    f76e:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    f772:	48 89 ce             	mov    %rcx,%rsi
    f775:	48 89 c7             	mov    %rax,%rdi
    f778:	e8 93 ab ff ff       	call   a310 <MPI_Comm_create@plt>
  return newcomm;
    f77d:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f781:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f785:	48 89 d6             	mov    %rdx,%rsi
    f788:	48 89 c7             	mov    %rax,%rdi
    f78b:	e8 22 c6 ff ff       	call   bdb2 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
}
    f790:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f794:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f79b:	00 00 
    f79d:	74 05                	je     f7a4 <_ZNK3MPI9Intercomm6CreateERKNS_5GroupE+0x74>
    f79f:	e8 3c ad ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f7a4:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f7a8:	c9                   	leave  
    f7a9:	c3                   	ret    

000000000000f7aa <_ZNK3MPI9Intercomm5SplitEii>:

inline MPI::Intercomm
MPI::Intercomm::Split(int color, int key) const
{
    f7aa:	f3 0f 1e fa          	endbr64 
    f7ae:	55                   	push   %rbp
    f7af:	48 89 e5             	mov    %rsp,%rbp
    f7b2:	48 83 ec 30          	sub    $0x30,%rsp
    f7b6:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f7ba:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f7be:	89 55 dc             	mov    %edx,-0x24(%rbp)
    f7c1:	89 4d d8             	mov    %ecx,-0x28(%rbp)
    f7c4:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f7cb:	00 00 
    f7cd:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f7d1:	31 c0                	xor    %eax,%eax
  MPI_Comm newcomm;
  (void) MPI_Comm_split(mpi_comm, color, key, &newcomm);
    f7d3:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f7d7:	48 8b 40 08          	mov    0x8(%rax),%rax
    f7db:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    f7df:	8b 55 d8             	mov    -0x28(%rbp),%edx
    f7e2:	8b 75 dc             	mov    -0x24(%rbp),%esi
    f7e5:	48 89 c7             	mov    %rax,%rdi
    f7e8:	e8 a3 b0 ff ff       	call   a890 <MPI_Comm_split@plt>
  return newcomm;
    f7ed:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f7f1:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f7f5:	48 89 d6             	mov    %rdx,%rsi
    f7f8:	48 89 c7             	mov    %rax,%rdi
    f7fb:	e8 b2 c5 ff ff       	call   bdb2 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
}
    f800:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f804:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f80b:	00 00 
    f80d:	74 05                	je     f814 <_ZNK3MPI9Intercomm5SplitEii+0x6a>
    f80f:	e8 cc ac ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f814:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f818:	c9                   	leave  
    f819:	c3                   	ret    

000000000000f81a <_ZNK3MPI5Group8Get_sizeEv>:
// Groups, Contexts, and Communicators
//

inline int
MPI::Group::Get_size() const
{
    f81a:	f3 0f 1e fa          	endbr64 
    f81e:	55                   	push   %rbp
    f81f:	48 89 e5             	mov    %rsp,%rbp
    f822:	48 83 ec 20          	sub    $0x20,%rsp
    f826:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f82a:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f831:	00 00 
    f833:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f837:	31 c0                	xor    %eax,%eax
  int size;
  (void)MPI_Group_size(mpi_group, &size);
    f839:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f83d:	48 8b 40 08          	mov    0x8(%rax),%rax
    f841:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    f845:	48 89 d6             	mov    %rdx,%rsi
    f848:	48 89 c7             	mov    %rax,%rdi
    f84b:	e8 80 a6 ff ff       	call   9ed0 <MPI_Group_size@plt>
  return size;
    f850:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    f853:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    f857:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f85e:	00 00 
    f860:	74 05                	je     f867 <_ZNK3MPI5Group8Get_sizeEv+0x4d>
    f862:	e8 79 ac ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f867:	c9                   	leave  
    f868:	c3                   	ret    
    f869:	90                   	nop

000000000000f86a <_ZNK3MPI5Group8Get_rankEv>:

inline int
MPI::Group::Get_rank() const
{
    f86a:	f3 0f 1e fa          	endbr64 
    f86e:	55                   	push   %rbp
    f86f:	48 89 e5             	mov    %rsp,%rbp
    f872:	48 83 ec 20          	sub    $0x20,%rsp
    f876:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f87a:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f881:	00 00 
    f883:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f887:	31 c0                	xor    %eax,%eax
  int myrank;
  (void)MPI_Group_rank(mpi_group, &myrank);
    f889:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f88d:	48 8b 40 08          	mov    0x8(%rax),%rax
    f891:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    f895:	48 89 d6             	mov    %rdx,%rsi
    f898:	48 89 c7             	mov    %rax,%rdi
    f89b:	e8 90 b1 ff ff       	call   aa30 <MPI_Group_rank@plt>
  return myrank;
    f8a0:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    f8a3:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    f8a7:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f8ae:	00 00 
    f8b0:	74 05                	je     f8b7 <_ZNK3MPI5Group8Get_rankEv+0x4d>
    f8b2:	e8 29 ac ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f8b7:	c9                   	leave  
    f8b8:	c3                   	ret    
    f8b9:	90                   	nop

000000000000f8ba <_ZNK3MPI5Group4InclEiPKi>:
  return newgroup;
}

inline MPI::Group
MPI::Group::Incl(int n, const int ranks[]) const
{
    f8ba:	f3 0f 1e fa          	endbr64 
    f8be:	55                   	push   %rbp
    f8bf:	48 89 e5             	mov    %rsp,%rbp
    f8c2:	48 83 ec 30          	sub    $0x30,%rsp
    f8c6:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f8ca:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f8ce:	89 55 dc             	mov    %edx,-0x24(%rbp)
    f8d1:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    f8d5:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f8dc:	00 00 
    f8de:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f8e2:	31 c0                	xor    %eax,%eax
  MPI_Group newgroup;
  (void)MPI_Group_incl(mpi_group, n, const_cast<int *>(ranks), &newgroup);
    f8e4:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f8e8:	48 8b 40 08          	mov    0x8(%rax),%rax
    f8ec:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    f8f0:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    f8f4:	8b 75 dc             	mov    -0x24(%rbp),%esi
    f8f7:	48 89 c7             	mov    %rax,%rdi
    f8fa:	e8 31 a9 ff ff       	call   a230 <MPI_Group_incl@plt>
  return newgroup;
    f8ff:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f903:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f907:	48 89 d6             	mov    %rdx,%rsi
    f90a:	48 89 c7             	mov    %rax,%rdi
    f90d:	e8 7a c1 ff ff       	call   ba8c <_ZN3MPI5GroupC1EP12ompi_group_t>
}
    f912:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f916:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f91d:	00 00 
    f91f:	74 05                	je     f926 <_ZNK3MPI5Group4InclEiPKi+0x6c>
    f921:	e8 ba ab ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f926:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f92a:	c9                   	leave  
    f92b:	c3                   	ret    

000000000000f92c <_ZNK3MPI5Group4ExclEiPKi>:

inline MPI::Group
MPI::Group::Excl(int n, const int ranks[]) const
{
    f92c:	f3 0f 1e fa          	endbr64 
    f930:	55                   	push   %rbp
    f931:	48 89 e5             	mov    %rsp,%rbp
    f934:	48 83 ec 30          	sub    $0x30,%rsp
    f938:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f93c:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f940:	89 55 dc             	mov    %edx,-0x24(%rbp)
    f943:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    f947:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f94e:	00 00 
    f950:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f954:	31 c0                	xor    %eax,%eax
  MPI_Group newgroup;
  (void)MPI_Group_excl(mpi_group, n, const_cast<int *>(ranks), &newgroup);
    f956:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f95a:	48 8b 40 08          	mov    0x8(%rax),%rax
    f95e:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    f962:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    f966:	8b 75 dc             	mov    -0x24(%rbp),%esi
    f969:	48 89 c7             	mov    %rax,%rdi
    f96c:	e8 8f a4 ff ff       	call   9e00 <MPI_Group_excl@plt>
  return newgroup;
    f971:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f975:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f979:	48 89 d6             	mov    %rdx,%rsi
    f97c:	48 89 c7             	mov    %rax,%rdi
    f97f:	e8 08 c1 ff ff       	call   ba8c <_ZN3MPI5GroupC1EP12ompi_group_t>
}
    f984:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f988:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f98f:	00 00 
    f991:	74 05                	je     f998 <_ZNK3MPI5Group4ExclEiPKi+0x6c>
    f993:	e8 48 ab ff ff       	call   a4e0 <__stack_chk_fail@plt>
    f998:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f99c:	c9                   	leave  
    f99d:	c3                   	ret    

000000000000f99e <_ZNK3MPI5Group10Range_inclEiPA3_Ki>:

inline MPI::Group
MPI::Group::Range_incl(int n, const int ranges[][3]) const
{
    f99e:	f3 0f 1e fa          	endbr64 
    f9a2:	55                   	push   %rbp
    f9a3:	48 89 e5             	mov    %rsp,%rbp
    f9a6:	48 83 ec 30          	sub    $0x30,%rsp
    f9aa:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f9ae:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f9b2:	89 55 dc             	mov    %edx,-0x24(%rbp)
    f9b5:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    f9b9:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f9c0:	00 00 
    f9c2:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f9c6:	31 c0                	xor    %eax,%eax
  MPI_Group newgroup;
  (void)MPI_Group_range_incl(mpi_group, n,
    f9c8:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f9cc:	48 8b 40 08          	mov    0x8(%rax),%rax
    f9d0:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    f9d4:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    f9d8:	8b 75 dc             	mov    -0x24(%rbp),%esi
    f9db:	48 89 c7             	mov    %rax,%rdi
    f9de:	e8 bd a4 ff ff       	call   9ea0 <MPI_Group_range_incl@plt>
                             const_cast<int(*)[3]>(ranges),
#else
                             (int(*)[3]) ranges,
#endif
                             &newgroup);
  return newgroup;
    f9e3:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f9e7:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f9eb:	48 89 d6             	mov    %rdx,%rsi
    f9ee:	48 89 c7             	mov    %rax,%rdi
    f9f1:	e8 96 c0 ff ff       	call   ba8c <_ZN3MPI5GroupC1EP12ompi_group_t>
}
    f9f6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f9fa:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    fa01:	00 00 
    fa03:	74 05                	je     fa0a <_ZNK3MPI5Group10Range_inclEiPA3_Ki+0x6c>
    fa05:	e8 d6 aa ff ff       	call   a4e0 <__stack_chk_fail@plt>
    fa0a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fa0e:	c9                   	leave  
    fa0f:	c3                   	ret    

000000000000fa10 <_ZNK3MPI5Group10Range_exclEiPA3_Ki>:

inline MPI::Group
MPI::Group::Range_excl(int n, const int ranges[][3]) const
{
    fa10:	f3 0f 1e fa          	endbr64 
    fa14:	55                   	push   %rbp
    fa15:	48 89 e5             	mov    %rsp,%rbp
    fa18:	48 83 ec 30          	sub    $0x30,%rsp
    fa1c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fa20:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    fa24:	89 55 dc             	mov    %edx,-0x24(%rbp)
    fa27:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    fa2b:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    fa32:	00 00 
    fa34:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    fa38:	31 c0                	xor    %eax,%eax
  MPI_Group newgroup;
  (void)MPI_Group_range_excl(mpi_group, n,
    fa3a:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    fa3e:	48 8b 40 08          	mov    0x8(%rax),%rax
    fa42:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    fa46:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    fa4a:	8b 75 dc             	mov    -0x24(%rbp),%esi
    fa4d:	48 89 c7             	mov    %rax,%rdi
    fa50:	e8 3b ab ff ff       	call   a590 <MPI_Group_range_excl@plt>
                             const_cast<int(*)[3]>(ranges),
#else
                             (int(*)[3]) ranges,
#endif
                             &newgroup);
  return newgroup;
    fa55:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    fa59:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fa5d:	48 89 d6             	mov    %rdx,%rsi
    fa60:	48 89 c7             	mov    %rax,%rdi
    fa63:	e8 24 c0 ff ff       	call   ba8c <_ZN3MPI5GroupC1EP12ompi_group_t>
}
    fa68:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fa6c:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    fa73:	00 00 
    fa75:	74 05                	je     fa7c <_ZNK3MPI5Group10Range_exclEiPA3_Ki+0x6c>
    fa77:	e8 64 aa ff ff       	call   a4e0 <__stack_chk_fail@plt>
    fa7c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fa80:	c9                   	leave  
    fa81:	c3                   	ret    

000000000000fa82 <_ZN3MPI5Group4FreeEv>:

inline void
MPI::Group::Free()
{
    fa82:	f3 0f 1e fa          	endbr64 
    fa86:	55                   	push   %rbp
    fa87:	48 89 e5             	mov    %rsp,%rbp
    fa8a:	48 83 ec 10          	sub    $0x10,%rsp
    fa8e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
  (void)MPI_Group_free(&mpi_group);
    fa92:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fa96:	48 83 c0 08          	add    $0x8,%rax
    fa9a:	48 89 c7             	mov    %rax,%rdi
    fa9d:	e8 ee ae ff ff       	call   a990 <MPI_Group_free@plt>
}
    faa2:	90                   	nop
    faa3:	c9                   	leave  
    faa4:	c3                   	ret    
    faa5:	90                   	nop

000000000000faa6 <_ZN3MPI2OpD1Ev>:
inline
MPI::Op::Op(const MPI::Op& op)
  : mpi_op(op.mpi_op) { }

inline
MPI::Op::~Op()
    faa6:	f3 0f 1e fa          	endbr64 
    faaa:	55                   	push   %rbp
    faab:	48 89 e5             	mov    %rsp,%rbp
    faae:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fab2:	48 8d 15 3f bc 00 00 	lea    0xbc3f(%rip),%rdx        # 1b6f8 <_ZTVN3MPI2OpE+0x10>
    fab9:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fabd:	48 89 10             	mov    %rdx,(%rax)
{
#if 0
    mpi_op = MPI_OP_NULL;
    op_user_function = 0;
#endif
}
    fac0:	90                   	nop
    fac1:	5d                   	pop    %rbp
    fac2:	c3                   	ret    
    fac3:	90                   	nop

000000000000fac4 <_ZN3MPI2OpD0Ev>:
MPI::Op::~Op()
    fac4:	f3 0f 1e fa          	endbr64 
    fac8:	55                   	push   %rbp
    fac9:	48 89 e5             	mov    %rsp,%rbp
    facc:	48 83 ec 10          	sub    $0x10,%rsp
    fad0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
}
    fad4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fad8:	48 89 c7             	mov    %rax,%rdi
    fadb:	e8 c6 ff ff ff       	call   faa6 <_ZN3MPI2OpD1Ev>
    fae0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fae4:	be 10 00 00 00       	mov    $0x10,%esi
    fae9:	48 89 c7             	mov    %rax,%rdi
    faec:	e8 0f a9 ff ff       	call   a400 <_ZdlPvm@plt>
    faf1:	c9                   	leave  
    faf2:	c3                   	ret    
    faf3:	90                   	nop

000000000000faf4 <_ZNK3MPI2OpcvP9ompi_op_tEv>:
// inter-language operability
inline MPI::Op&
MPI::Op::operator= (const MPI_Op &i) { mpi_op = i; return *this; }

inline
MPI::Op::operator MPI_Op () const { return mpi_op; }
    faf4:	f3 0f 1e fa          	endbr64 
    faf8:	55                   	push   %rbp
    faf9:	48 89 e5             	mov    %rsp,%rbp
    fafc:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fb00:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fb04:	48 8b 40 08          	mov    0x8(%rax),%rax
    fb08:	5d                   	pop    %rbp
    fb09:	c3                   	ret    

000000000000fb0a <_ZN3MPI2Op4InitEPFvPKvPviRKNS_8DatatypeEEb>:
// There is a lengthy comment in ompi/mpi/cxx/intercepts.cc explaining
// what this function is doing.  Please read it before modifying this
// function.
inline void
MPI::Op::Init(MPI::User_function *func, bool commute)
{
    fb0a:	f3 0f 1e fa          	endbr64 
    fb0e:	55                   	push   %rbp
    fb0f:	48 89 e5             	mov    %rsp,%rbp
    fb12:	48 83 ec 20          	sub    $0x20,%rsp
    fb16:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fb1a:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    fb1e:	89 d0                	mov    %edx,%eax
    fb20:	88 45 ec             	mov    %al,-0x14(%rbp)
    (void)MPI_Op_create((MPI_User_function*) ompi_mpi_cxx_op_intercept,
    fb23:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fb27:	48 8d 50 08          	lea    0x8(%rax),%rdx
    fb2b:	0f b6 45 ec          	movzbl -0x14(%rbp),%eax
    fb2f:	89 c6                	mov    %eax,%esi
    fb31:	48 8b 05 88 d4 00 00 	mov    0xd488(%rip),%rax        # 1cfc0 <ompi_mpi_cxx_op_intercept@Base>
    fb38:	48 89 c7             	mov    %rax,%rdi
    fb3b:	e8 80 a3 ff ff       	call   9ec0 <MPI_Op_create@plt>
                        (int) commute, &mpi_op);
    ompi_op_set_cxx_callback(mpi_op, (MPI_User_function*) func);
    fb40:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fb44:	48 8b 40 08          	mov    0x8(%rax),%rax
    fb48:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    fb4c:	48 89 d6             	mov    %rdx,%rsi
    fb4f:	48 89 c7             	mov    %rax,%rdi
    fb52:	e8 09 a9 ff ff       	call   a460 <ompi_op_set_cxx_callback@plt>
}
    fb57:	90                   	nop
    fb58:	c9                   	leave  
    fb59:	c3                   	ret    

000000000000fb5a <_ZN3MPI2Op4FreeEv>:


inline void
MPI::Op::Free()
{
    fb5a:	f3 0f 1e fa          	endbr64 
    fb5e:	55                   	push   %rbp
    fb5f:	48 89 e5             	mov    %rsp,%rbp
    fb62:	48 83 ec 10          	sub    $0x10,%rsp
    fb66:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    (void)MPI_Op_free(&mpi_op);
    fb6a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fb6e:	48 83 c0 08          	add    $0x8,%rax
    fb72:	48 89 c7             	mov    %rax,%rdi
    fb75:	e8 a6 ac ff ff       	call   a820 <MPI_Op_free@plt>
}
    fb7a:	90                   	nop
    fb7b:	c9                   	leave  
    fb7c:	c3                   	ret    
    fb7d:	90                   	nop

000000000000fb7e <_ZNK3MPI2Op12Reduce_localEPKvPviRKNS_8DatatypeE>:


inline void
MPI::Op::Reduce_local(const void *inbuf, void *inoutbuf, int count,
                      const MPI::Datatype& datatype) const
{
    fb7e:	f3 0f 1e fa          	endbr64 
    fb82:	55                   	push   %rbp
    fb83:	48 89 e5             	mov    %rsp,%rbp
    fb86:	53                   	push   %rbx
    fb87:	48 83 ec 38          	sub    $0x38,%rsp
    fb8b:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fb8f:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    fb93:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    fb97:	89 4d d4             	mov    %ecx,-0x2c(%rbp)
    fb9a:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    (void)MPI_Reduce_local(const_cast<void*>(inbuf), inoutbuf, count,
                           datatype, mpi_op);
    fb9e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fba2:	48 8b 58 08          	mov    0x8(%rax),%rbx
    (void)MPI_Reduce_local(const_cast<void*>(inbuf), inoutbuf, count,
    fba6:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    fbaa:	48 89 c7             	mov    %rax,%rdi
    fbad:	e8 cc bc ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    fbb2:	48 89 c1             	mov    %rax,%rcx
    fbb5:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    fbb8:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    fbbc:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    fbc0:	49 89 d8             	mov    %rbx,%r8
    fbc3:	48 89 c7             	mov    %rax,%rdi
    fbc6:	e8 55 a9 ff ff       	call   a520 <MPI_Reduce_local@plt>
}
    fbcb:	90                   	nop
    fbcc:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    fbd0:	c9                   	leave  
    fbd1:	c3                   	ret    

000000000000fbd2 <_ZNK3MPI2Op14Is_commutativeEv>:


inline bool
MPI::Op::Is_commutative(void) const
{
    fbd2:	f3 0f 1e fa          	endbr64 
    fbd6:	55                   	push   %rbp
    fbd7:	48 89 e5             	mov    %rsp,%rbp
    fbda:	48 83 ec 20          	sub    $0x20,%rsp
    fbde:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fbe2:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    fbe9:	00 00 
    fbeb:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    fbef:	31 c0                	xor    %eax,%eax
    int commute;
    (void)MPI_Op_commutative(mpi_op, &commute);
    fbf1:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fbf5:	48 8b 40 08          	mov    0x8(%rax),%rax
    fbf9:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    fbfd:	48 89 d6             	mov    %rdx,%rsi
    fc00:	48 89 c7             	mov    %rax,%rdi
    fc03:	e8 48 a7 ff ff       	call   a350 <MPI_Op_commutative@plt>
    return (bool) commute;
    fc08:	8b 45 f4             	mov    -0xc(%rbp),%eax
    fc0b:	85 c0                	test   %eax,%eax
    fc0d:	0f 95 c0             	setne  %al
}
    fc10:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    fc14:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    fc1b:	00 00 
    fc1d:	74 05                	je     fc24 <_ZNK3MPI2Op14Is_commutativeEv+0x52>
    fc1f:	e8 bc a8 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    fc24:	c9                   	leave  
    fc25:	c3                   	ret    

000000000000fc26 <_ZN3MPI10Errhandler4FreeEv>:

#endif

inline void
MPI::Errhandler::Free()
{
    fc26:	f3 0f 1e fa          	endbr64 
    fc2a:	55                   	push   %rbp
    fc2b:	48 89 e5             	mov    %rsp,%rbp
    fc2e:	48 83 ec 10          	sub    $0x10,%rsp
    fc32:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
  (void)MPI_Errhandler_free(&mpi_errhandler);
    fc36:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fc3a:	48 83 c0 08          	add    $0x8,%rax
    fc3e:	48 89 c7             	mov    %rax,%rdi
    fc41:	e8 aa a7 ff ff       	call   a3f0 <MPI_Errhandler_free@plt>
}
    fc46:	90                   	nop
    fc47:	c9                   	leave  
    fc48:	c3                   	ret    
    fc49:	90                   	nop

000000000000fc4a <_ZNK3MPI6Status9Get_countERKNS_8DatatypeE>:
// Point-to-Point Communication
//

inline int
MPI::Status::Get_count(const MPI::Datatype& datatype) const
{
    fc4a:	f3 0f 1e fa          	endbr64 
    fc4e:	55                   	push   %rbp
    fc4f:	48 89 e5             	mov    %rsp,%rbp
    fc52:	48 83 ec 20          	sub    $0x20,%rsp
    fc56:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fc5a:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    fc5e:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    fc65:	00 00 
    fc67:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    fc6b:	31 c0                	xor    %eax,%eax
  int count;
  (void)MPI_Get_count(const_cast<MPI_Status*>(&mpi_status), datatype, &count);
    fc6d:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    fc71:	48 89 c7             	mov    %rax,%rdi
    fc74:	e8 05 bc ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    fc79:	48 89 c6             	mov    %rax,%rsi
    fc7c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fc80:	48 8d 48 08          	lea    0x8(%rax),%rcx
    fc84:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
    fc88:	48 89 c2             	mov    %rax,%rdx
    fc8b:	48 89 cf             	mov    %rcx,%rdi
    fc8e:	e8 bd aa ff ff       	call   a750 <MPI_Get_count@plt>
  return count;
    fc93:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    fc96:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    fc9a:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    fca1:	00 00 
    fca3:	74 05                	je     fcaa <_ZNK3MPI6Status9Get_countERKNS_8DatatypeE+0x60>
    fca5:	e8 36 a8 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    fcaa:	c9                   	leave  
    fcab:	c3                   	ret    

000000000000fcac <_ZNK3MPI6Status12Is_cancelledEv>:

inline bool
MPI::Status::Is_cancelled() const
{
    fcac:	f3 0f 1e fa          	endbr64 
    fcb0:	55                   	push   %rbp
    fcb1:	48 89 e5             	mov    %rsp,%rbp
    fcb4:	48 83 ec 20          	sub    $0x20,%rsp
    fcb8:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fcbc:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    fcc3:	00 00 
    fcc5:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    fcc9:	31 c0                	xor    %eax,%eax
  int t;
  (void)MPI_Test_cancelled(const_cast<MPI_Status*>(&mpi_status), &t);
    fccb:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fccf:	48 8d 50 08          	lea    0x8(%rax),%rdx
    fcd3:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
    fcd7:	48 89 c6             	mov    %rax,%rsi
    fcda:	48 89 d7             	mov    %rdx,%rdi
    fcdd:	e8 9e a6 ff ff       	call   a380 <MPI_Test_cancelled@plt>
  return OPAL_INT_TO_BOOL(t);
    fce2:	8b 45 f4             	mov    -0xc(%rbp),%eax
    fce5:	85 c0                	test   %eax,%eax
    fce7:	0f 95 c0             	setne  %al
}
    fcea:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    fcee:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    fcf5:	00 00 
    fcf7:	74 05                	je     fcfe <_ZNK3MPI6Status12Is_cancelledEv+0x52>
    fcf9:	e8 e2 a7 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    fcfe:	c9                   	leave  
    fcff:	c3                   	ret    

000000000000fd00 <_ZNK3MPI6Status12Get_elementsERKNS_8DatatypeE>:

inline int
MPI::Status::Get_elements(const MPI::Datatype& datatype) const
{
    fd00:	f3 0f 1e fa          	endbr64 
    fd04:	55                   	push   %rbp
    fd05:	48 89 e5             	mov    %rsp,%rbp
    fd08:	48 83 ec 20          	sub    $0x20,%rsp
    fd0c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fd10:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    fd14:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    fd1b:	00 00 
    fd1d:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    fd21:	31 c0                	xor    %eax,%eax
  int count;
  (void)MPI_Get_elements(const_cast<MPI_Status*>(&mpi_status), datatype, &count);
    fd23:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    fd27:	48 89 c7             	mov    %rax,%rdi
    fd2a:	e8 4f bb ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    fd2f:	48 89 c6             	mov    %rax,%rsi
    fd32:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fd36:	48 8d 48 08          	lea    0x8(%rax),%rcx
    fd3a:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
    fd3e:	48 89 c2             	mov    %rax,%rdx
    fd41:	48 89 cf             	mov    %rcx,%rdi
    fd44:	e8 37 a7 ff ff       	call   a480 <MPI_Get_elements@plt>
  return count;
    fd49:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    fd4c:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    fd50:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    fd57:	00 00 
    fd59:	74 05                	je     fd60 <_ZNK3MPI6Status12Get_elementsERKNS_8DatatypeE+0x60>
    fd5b:	e8 80 a7 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    fd60:	c9                   	leave  
    fd61:	c3                   	ret    

000000000000fd62 <_ZNK3MPI6Status10Get_sourceEv>:
//
// Status Access
//
inline int
MPI::Status::Get_source() const
{
    fd62:	f3 0f 1e fa          	endbr64 
    fd66:	55                   	push   %rbp
    fd67:	48 89 e5             	mov    %rsp,%rbp
    fd6a:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
  int source;
  source = mpi_status.MPI_SOURCE;
    fd6e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fd72:	8b 40 08             	mov    0x8(%rax),%eax
    fd75:	89 45 fc             	mov    %eax,-0x4(%rbp)
  return source;
    fd78:	8b 45 fc             	mov    -0x4(%rbp),%eax
}
    fd7b:	5d                   	pop    %rbp
    fd7c:	c3                   	ret    
    fd7d:	90                   	nop

000000000000fd7e <_ZN3MPI6Status10Set_sourceEi>:

inline void
MPI::Status::Set_source(int source)
{
    fd7e:	f3 0f 1e fa          	endbr64 
    fd82:	55                   	push   %rbp
    fd83:	48 89 e5             	mov    %rsp,%rbp
    fd86:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fd8a:	89 75 f4             	mov    %esi,-0xc(%rbp)
  mpi_status.MPI_SOURCE = source;
    fd8d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fd91:	8b 55 f4             	mov    -0xc(%rbp),%edx
    fd94:	89 50 08             	mov    %edx,0x8(%rax)
}
    fd97:	90                   	nop
    fd98:	5d                   	pop    %rbp
    fd99:	c3                   	ret    

000000000000fd9a <_ZNK3MPI6Status7Get_tagEv>:

inline int
MPI::Status::Get_tag() const
{
    fd9a:	f3 0f 1e fa          	endbr64 
    fd9e:	55                   	push   %rbp
    fd9f:	48 89 e5             	mov    %rsp,%rbp
    fda2:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
  int tag;
  tag = mpi_status.MPI_TAG;
    fda6:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fdaa:	8b 40 0c             	mov    0xc(%rax),%eax
    fdad:	89 45 fc             	mov    %eax,-0x4(%rbp)
  return tag;
    fdb0:	8b 45 fc             	mov    -0x4(%rbp),%eax
}
    fdb3:	5d                   	pop    %rbp
    fdb4:	c3                   	ret    
    fdb5:	90                   	nop

000000000000fdb6 <_ZN3MPI6Status7Set_tagEi>:

inline void
MPI::Status::Set_tag(int tag)
{
    fdb6:	f3 0f 1e fa          	endbr64 
    fdba:	55                   	push   %rbp
    fdbb:	48 89 e5             	mov    %rsp,%rbp
    fdbe:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fdc2:	89 75 f4             	mov    %esi,-0xc(%rbp)
  mpi_status.MPI_TAG = tag;
    fdc5:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fdc9:	8b 55 f4             	mov    -0xc(%rbp),%edx
    fdcc:	89 50 0c             	mov    %edx,0xc(%rax)
}
    fdcf:	90                   	nop
    fdd0:	5d                   	pop    %rbp
    fdd1:	c3                   	ret    

000000000000fdd2 <_ZNK3MPI6Status9Get_errorEv>:

inline int
MPI::Status::Get_error() const
{
    fdd2:	f3 0f 1e fa          	endbr64 
    fdd6:	55                   	push   %rbp
    fdd7:	48 89 e5             	mov    %rsp,%rbp
    fdda:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
  int error;
  error = mpi_status.MPI_ERROR;
    fdde:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fde2:	8b 40 10             	mov    0x10(%rax),%eax
    fde5:	89 45 fc             	mov    %eax,-0x4(%rbp)
  return error;
    fde8:	8b 45 fc             	mov    -0x4(%rbp),%eax
}
    fdeb:	5d                   	pop    %rbp
    fdec:	c3                   	ret    
    fded:	90                   	nop

000000000000fdee <_ZN3MPI6Status9Set_errorEi>:

inline void
MPI::Status::Set_error(int error)
{
    fdee:	f3 0f 1e fa          	endbr64 
    fdf2:	55                   	push   %rbp
    fdf3:	48 89 e5             	mov    %rsp,%rbp
    fdf6:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fdfa:	89 75 f4             	mov    %esi,-0xc(%rbp)
  mpi_status.MPI_ERROR = error;
    fdfd:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fe01:	8b 55 f4             	mov    -0xc(%rbp),%edx
    fe04:	89 50 10             	mov    %edx,0x10(%rax)
}
    fe07:	90                   	nop
    fe08:	5d                   	pop    %rbp
    fe09:	c3                   	ret    

000000000000fe0a <_ZN3MPI6Status12Set_elementsERKNS_8DatatypeEi>:

inline void
MPI::Status::Set_elements(const MPI::Datatype& datatype, int count)
{
    fe0a:	f3 0f 1e fa          	endbr64 
    fe0e:	55                   	push   %rbp
    fe0f:	48 89 e5             	mov    %rsp,%rbp
    fe12:	48 83 ec 20          	sub    $0x20,%rsp
    fe16:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fe1a:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    fe1e:	89 55 ec             	mov    %edx,-0x14(%rbp)
    MPI_Status_set_elements(&mpi_status, datatype, count);
    fe21:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    fe25:	48 89 c7             	mov    %rax,%rdi
    fe28:	e8 51 ba ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    fe2d:	48 89 c6             	mov    %rax,%rsi
    fe30:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fe34:	48 8d 48 08          	lea    0x8(%rax),%rcx
    fe38:	8b 45 ec             	mov    -0x14(%rbp),%eax
    fe3b:	89 c2                	mov    %eax,%edx
    fe3d:	48 89 cf             	mov    %rcx,%rdi
    fe40:	e8 ab a7 ff ff       	call   a5f0 <MPI_Status_set_elements@plt>
}
    fe45:	90                   	nop
    fe46:	c9                   	leave  
    fe47:	c3                   	ret    

000000000000fe48 <_ZN3MPI6Status13Set_cancelledEb>:

inline void
MPI::Status::Set_cancelled(bool flag)
{
    fe48:	f3 0f 1e fa          	endbr64 
    fe4c:	55                   	push   %rbp
    fe4d:	48 89 e5             	mov    %rsp,%rbp
    fe50:	48 83 ec 10          	sub    $0x10,%rsp
    fe54:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fe58:	89 f0                	mov    %esi,%eax
    fe5a:	88 45 f4             	mov    %al,-0xc(%rbp)
    MPI_Status_set_cancelled(&mpi_status, (int) flag);
    fe5d:	0f b6 45 f4          	movzbl -0xc(%rbp),%eax
    fe61:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    fe65:	48 83 c2 08          	add    $0x8,%rdx
    fe69:	89 c6                	mov    %eax,%esi
    fe6b:	48 89 d7             	mov    %rdx,%rdi
    fe6e:	e8 2d ac ff ff       	call   aaa0 <MPI_Status_set_cancelled@plt>
}
    fe73:	90                   	nop
    fe74:	c9                   	leave  
    fe75:	c3                   	ret    

000000000000fe76 <_ZN3MPI4Info6DeleteEPKc>:
  return newinfo;
}

inline void
MPI::Info::Delete(const char* key)
{
    fe76:	f3 0f 1e fa          	endbr64 
    fe7a:	55                   	push   %rbp
    fe7b:	48 89 e5             	mov    %rsp,%rbp
    fe7e:	48 83 ec 10          	sub    $0x10,%rsp
    fe82:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fe86:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   (void)MPI_Info_delete(mpi_info, const_cast<char *>(key));
    fe8a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fe8e:	48 8b 40 08          	mov    0x8(%rax),%rax
    fe92:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    fe96:	48 89 d6             	mov    %rdx,%rsi
    fe99:	48 89 c7             	mov    %rax,%rdi
    fe9c:	e8 8f a8 ff ff       	call   a730 <MPI_Info_delete@plt>
}
    fea1:	90                   	nop
    fea2:	c9                   	leave  
    fea3:	c3                   	ret    

000000000000fea4 <_ZNK3MPI4Info3DupEv>:

inline MPI::Info
MPI::Info::Dup() const
{
    fea4:	f3 0f 1e fa          	endbr64 
    fea8:	55                   	push   %rbp
    fea9:	48 89 e5             	mov    %rsp,%rbp
    feac:	48 83 ec 20          	sub    $0x20,%rsp
    feb0:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    feb4:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    feb8:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    febf:	00 00 
    fec1:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    fec5:	31 c0                	xor    %eax,%eax
  MPI_Info newinfo;
  (void)MPI_Info_dup(mpi_info, &newinfo);
    fec7:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    fecb:	48 8b 40 08          	mov    0x8(%rax),%rax
    fecf:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    fed3:	48 89 d6             	mov    %rdx,%rsi
    fed6:	48 89 c7             	mov    %rax,%rdi
    fed9:	e8 72 9f ff ff       	call   9e50 <MPI_Info_dup@plt>
  return newinfo;
    fede:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    fee2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fee6:	48 89 d6             	mov    %rdx,%rsi
    fee9:	48 89 c7             	mov    %rax,%rdi
    feec:	e8 f9 be ff ff       	call   bdea <_ZN3MPI4InfoC1EP11ompi_info_t>
}
    fef1:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fef5:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    fefc:	00 00 
    fefe:	74 05                	je     ff05 <_ZNK3MPI4Info3DupEv+0x61>
    ff00:	e8 db a5 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    ff05:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ff09:	c9                   	leave  
    ff0a:	c3                   	ret    
    ff0b:	90                   	nop

000000000000ff0c <_ZN3MPI4Info4FreeEv>:

inline void
MPI::Info::Free()
{
    ff0c:	f3 0f 1e fa          	endbr64 
    ff10:	55                   	push   %rbp
    ff11:	48 89 e5             	mov    %rsp,%rbp
    ff14:	48 83 ec 10          	sub    $0x10,%rsp
    ff18:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
  (void) MPI_Info_free(&mpi_info);
    ff1c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ff20:	48 83 c0 08          	add    $0x8,%rax
    ff24:	48 89 c7             	mov    %rax,%rdi
    ff27:	e8 94 ab ff ff       	call   aac0 <MPI_Info_free@plt>
}
    ff2c:	90                   	nop
    ff2d:	c9                   	leave  
    ff2e:	c3                   	ret    
    ff2f:	90                   	nop

000000000000ff30 <_ZNK3MPI4Info3GetEPKciPc>:

inline bool
MPI::Info::Get(const char* key, int valuelen, char* value) const
{
    ff30:	f3 0f 1e fa          	endbr64 
    ff34:	55                   	push   %rbp
    ff35:	48 89 e5             	mov    %rsp,%rbp
    ff38:	48 83 ec 30          	sub    $0x30,%rsp
    ff3c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    ff40:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    ff44:	89 55 dc             	mov    %edx,-0x24(%rbp)
    ff47:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    ff4b:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ff52:	00 00 
    ff54:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    ff58:	31 c0                	xor    %eax,%eax
  int flag;
  (void)MPI_Info_get(mpi_info, const_cast<char *>(key), valuelen, value, &flag);
    ff5a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ff5e:	48 8b 40 08          	mov    0x8(%rax),%rax
    ff62:	48 8d 7d f4          	lea    -0xc(%rbp),%rdi
    ff66:	48 8b 4d d0          	mov    -0x30(%rbp),%rcx
    ff6a:	8b 55 dc             	mov    -0x24(%rbp),%edx
    ff6d:	48 8b 75 e0          	mov    -0x20(%rbp),%rsi
    ff71:	49 89 f8             	mov    %rdi,%r8
    ff74:	48 89 c7             	mov    %rax,%rdi
    ff77:	e8 54 a1 ff ff       	call   a0d0 <MPI_Info_get@plt>
  return OPAL_INT_TO_BOOL(flag);
    ff7c:	8b 45 f4             	mov    -0xc(%rbp),%eax
    ff7f:	85 c0                	test   %eax,%eax
    ff81:	0f 95 c0             	setne  %al
}
    ff84:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    ff88:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    ff8f:	00 00 
    ff91:	74 05                	je     ff98 <_ZNK3MPI4Info3GetEPKciPc+0x68>
    ff93:	e8 48 a5 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    ff98:	c9                   	leave  
    ff99:	c3                   	ret    

000000000000ff9a <_ZNK3MPI4Info9Get_nkeysEv>:

inline int
MPI::Info::Get_nkeys() const
{
    ff9a:	f3 0f 1e fa          	endbr64 
    ff9e:	55                   	push   %rbp
    ff9f:	48 89 e5             	mov    %rsp,%rbp
    ffa2:	48 83 ec 20          	sub    $0x20,%rsp
    ffa6:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    ffaa:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ffb1:	00 00 
    ffb3:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    ffb7:	31 c0                	xor    %eax,%eax
  int nkeys;
  MPI_Info_get_nkeys(mpi_info, &nkeys);
    ffb9:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ffbd:	48 8b 40 08          	mov    0x8(%rax),%rax
    ffc1:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    ffc5:	48 89 d6             	mov    %rdx,%rsi
    ffc8:	48 89 c7             	mov    %rax,%rdi
    ffcb:	e8 e0 a4 ff ff       	call   a4b0 <MPI_Info_get_nkeys@plt>
  return nkeys;
    ffd0:	8b 45 f4             	mov    -0xc(%rbp),%eax
}
    ffd3:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    ffd7:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    ffde:	00 00 
    ffe0:	74 05                	je     ffe7 <_ZNK3MPI4Info9Get_nkeysEv+0x4d>
    ffe2:	e8 f9 a4 ff ff       	call   a4e0 <__stack_chk_fail@plt>
    ffe7:	c9                   	leave  
    ffe8:	c3                   	ret    
    ffe9:	90                   	nop

000000000000ffea <_ZNK3MPI4Info10Get_nthkeyEiPc>:

inline void
MPI::Info::Get_nthkey(int n, char* key) const
{
    ffea:	f3 0f 1e fa          	endbr64 
    ffee:	55                   	push   %rbp
    ffef:	48 89 e5             	mov    %rsp,%rbp
    fff2:	48 83 ec 20          	sub    $0x20,%rsp
    fff6:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fffa:	89 75 f4             	mov    %esi,-0xc(%rbp)
    fffd:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
  (void) MPI_Info_get_nthkey(mpi_info, n, key);
   10001:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10005:	48 8b 40 08          	mov    0x8(%rax),%rax
   10009:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   1000d:	8b 4d f4             	mov    -0xc(%rbp),%ecx
   10010:	89 ce                	mov    %ecx,%esi
   10012:	48 89 c7             	mov    %rax,%rdi
   10015:	e8 56 a7 ff ff       	call   a770 <MPI_Info_get_nthkey@plt>
}
   1001a:	90                   	nop
   1001b:	c9                   	leave  
   1001c:	c3                   	ret    
   1001d:	90                   	nop

000000000001001e <_ZNK3MPI4Info12Get_valuelenEPKcRi>:

inline bool
MPI::Info::Get_valuelen(const char* key, int& valuelen) const
{
   1001e:	f3 0f 1e fa          	endbr64 
   10022:	55                   	push   %rbp
   10023:	48 89 e5             	mov    %rsp,%rbp
   10026:	48 83 ec 30          	sub    $0x30,%rsp
   1002a:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   1002e:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   10032:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
   10036:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   1003d:	00 00 
   1003f:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   10043:	31 c0                	xor    %eax,%eax
  int flag;
  (void) MPI_Info_get_valuelen(mpi_info, const_cast<char *>(key), &valuelen, &flag);
   10045:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10049:	48 8b 40 08          	mov    0x8(%rax),%rax
   1004d:	48 8d 4d f4          	lea    -0xc(%rbp),%rcx
   10051:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
   10055:	48 8b 75 e0          	mov    -0x20(%rbp),%rsi
   10059:	48 89 c7             	mov    %rax,%rdi
   1005c:	e8 8f a8 ff ff       	call   a8f0 <MPI_Info_get_valuelen@plt>
  return OPAL_INT_TO_BOOL(flag);
   10061:	8b 45 f4             	mov    -0xc(%rbp),%eax
   10064:	85 c0                	test   %eax,%eax
   10066:	0f 95 c0             	setne  %al
}
   10069:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   1006d:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
   10074:	00 00 
   10076:	74 05                	je     1007d <_ZNK3MPI4Info12Get_valuelenEPKcRi+0x5f>
   10078:	e8 63 a4 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   1007d:	c9                   	leave  
   1007e:	c3                   	ret    
   1007f:	90                   	nop

0000000000010080 <_ZN3MPI4Info3SetEPKcS2_>:

inline void
MPI::Info::Set(const char* key, const char* value)
{
   10080:	f3 0f 1e fa          	endbr64 
   10084:	55                   	push   %rbp
   10085:	48 89 e5             	mov    %rsp,%rbp
   10088:	48 83 ec 20          	sub    $0x20,%rsp
   1008c:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10090:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   10094:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
  (void) MPI_Info_set(mpi_info, const_cast<char *>(key), const_cast<char *>(value));
   10098:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1009c:	48 8b 40 08          	mov    0x8(%rax),%rax
   100a0:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   100a4:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
   100a8:	48 89 ce             	mov    %rcx,%rsi
   100ab:	48 89 c7             	mov    %rax,%rdi
   100ae:	e8 1d 9f ff ff       	call   9fd0 <MPI_Info_set@plt>
}
   100b3:	90                   	nop
   100b4:	c9                   	leave  
   100b5:	c3                   	ret    

00000000000100b6 <_ZNK3MPI3Win14Get_errhandlerEv>:
//


inline MPI::Errhandler
MPI::Win:: Get_errhandler() const
{
   100b6:	f3 0f 1e fa          	endbr64 
   100ba:	55                   	push   %rbp
   100bb:	48 89 e5             	mov    %rsp,%rbp
   100be:	48 83 ec 20          	sub    $0x20,%rsp
   100c2:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   100c6:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   100ca:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   100d1:	00 00 
   100d3:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   100d7:	31 c0                	xor    %eax,%eax
    MPI_Errhandler errhandler;
    MPI_Win_get_errhandler(mpi_win, &errhandler);
   100d9:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   100dd:	48 8b 40 08          	mov    0x8(%rax),%rax
   100e1:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
   100e5:	48 89 d6             	mov    %rdx,%rsi
   100e8:	48 89 c7             	mov    %rax,%rdi
   100eb:	e8 70 a9 ff ff       	call   aa60 <MPI_Win_get_errhandler@plt>
    return errhandler;
   100f0:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   100f4:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   100f8:	48 89 d6             	mov    %rdx,%rsi
   100fb:	48 89 c7             	mov    %rax,%rdi
   100fe:	e8 81 bb ff ff       	call   bc84 <_ZN3MPI10ErrhandlerC1EP17ompi_errhandler_t>
}
   10103:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10107:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   1010e:	00 00 
   10110:	74 05                	je     10117 <_ZNK3MPI3Win14Get_errhandlerEv+0x61>
   10112:	e8 c9 a3 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   10117:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1011b:	c9                   	leave  
   1011c:	c3                   	ret    
   1011d:	90                   	nop

000000000001011e <_ZNK3MPI3Win14Set_errhandlerERKNS_10ErrhandlerE>:


inline void
MPI::Win::Set_errhandler(const MPI::Errhandler& errhandler) const
{
   1011e:	f3 0f 1e fa          	endbr64 
   10122:	55                   	push   %rbp
   10123:	48 89 e5             	mov    %rsp,%rbp
   10126:	48 83 ec 10          	sub    $0x10,%rsp
   1012a:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   1012e:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    (void)MPI_Win_set_errhandler(mpi_win, errhandler);
   10132:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   10136:	48 89 c7             	mov    %rax,%rdi
   10139:	e8 74 bb ff ff       	call   bcb2 <_ZNK3MPI10ErrhandlercvP17ompi_errhandler_tEv>
   1013e:	48 89 c2             	mov    %rax,%rdx
   10141:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10145:	48 8b 40 08          	mov    0x8(%rax),%rax
   10149:	48 89 d6             	mov    %rdx,%rsi
   1014c:	48 89 c7             	mov    %rax,%rdi
   1014f:	e8 3c a6 ff ff       	call   a790 <MPI_Win_set_errhandler@plt>
}
   10154:	90                   	nop
   10155:	c9                   	leave  
   10156:	c3                   	ret    
   10157:	90                   	nop

0000000000010158 <_ZNK3MPI3Win10AccumulateEPKviRKNS_8DatatypeEiliS5_RKNS_2OpE>:
MPI::Win::Accumulate(const void* origin_addr, int origin_count,
	       	const MPI::Datatype& origin_datatype, int target_rank,
	       	MPI::Aint target_disp, int target_count,
	       	const MPI::Datatype& target_datatype,
	       	const MPI::Op& op) const
{
   10158:	f3 0f 1e fa          	endbr64 
   1015c:	55                   	push   %rbp
   1015d:	48 89 e5             	mov    %rsp,%rbp
   10160:	41 55                	push   %r13
   10162:	41 54                	push   %r12
   10164:	53                   	push   %rbx
   10165:	48 83 ec 38          	sub    $0x38,%rsp
   10169:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   1016d:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
   10171:	89 55 cc             	mov    %edx,-0x34(%rbp)
   10174:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
   10178:	44 89 45 c8          	mov    %r8d,-0x38(%rbp)
   1017c:	4c 89 4d b8          	mov    %r9,-0x48(%rbp)
  (void) MPI_Accumulate(const_cast<void *>(origin_addr), origin_count, origin_datatype,
			target_rank, target_disp, target_count,
			target_datatype, op, mpi_win);
   10180:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10184:	4c 8b 68 08          	mov    0x8(%rax),%r13
  (void) MPI_Accumulate(const_cast<void *>(origin_addr), origin_count, origin_datatype,
   10188:	48 8b 45 20          	mov    0x20(%rbp),%rax
   1018c:	48 89 c7             	mov    %rax,%rdi
   1018f:	e8 60 f9 ff ff       	call   faf4 <_ZNK3MPI2OpcvP9ompi_op_tEv>
   10194:	49 89 c4             	mov    %rax,%r12
   10197:	48 8b 45 18          	mov    0x18(%rbp),%rax
   1019b:	48 89 c7             	mov    %rax,%rdi
   1019e:	e8 db b6 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
   101a3:	48 89 c3             	mov    %rax,%rbx
   101a6:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   101aa:	48 89 c7             	mov    %rax,%rdi
   101ad:	e8 cc b6 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
   101b2:	48 89 c7             	mov    %rax,%rdi
   101b5:	48 8b 4d b8          	mov    -0x48(%rbp),%rcx
   101b9:	8b 55 c8             	mov    -0x38(%rbp),%edx
   101bc:	8b 75 cc             	mov    -0x34(%rbp),%esi
   101bf:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   101c3:	48 83 ec 08          	sub    $0x8,%rsp
   101c7:	41 55                	push   %r13
   101c9:	41 54                	push   %r12
   101cb:	53                   	push   %rbx
   101cc:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
   101d0:	49 89 c8             	mov    %rcx,%r8
   101d3:	89 d1                	mov    %edx,%ecx
   101d5:	48 89 fa             	mov    %rdi,%rdx
   101d8:	48 89 c7             	mov    %rax,%rdi
   101db:	e8 90 a8 ff ff       	call   aa70 <MPI_Accumulate@plt>
   101e0:	48 83 c4 20          	add    $0x20,%rsp

}
   101e4:	90                   	nop
   101e5:	48 8d 65 e8          	lea    -0x18(%rbp),%rsp
   101e9:	5b                   	pop    %rbx
   101ea:	41 5c                	pop    %r12
   101ec:	41 5d                	pop    %r13
   101ee:	5d                   	pop    %rbp
   101ef:	c3                   	ret    

00000000000101f0 <_ZNK3MPI3Win8CompleteEv>:


inline void
MPI::Win::Complete() const
{
   101f0:	f3 0f 1e fa          	endbr64 
   101f4:	55                   	push   %rbp
   101f5:	48 89 e5             	mov    %rsp,%rbp
   101f8:	48 83 ec 10          	sub    $0x10,%rsp
   101fc:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
  (void) MPI_Win_complete(mpi_win);
   10200:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10204:	48 8b 40 08          	mov    0x8(%rax),%rax
   10208:	48 89 c7             	mov    %rax,%rdi
   1020b:	e8 60 a7 ff ff       	call   a970 <MPI_Win_complete@plt>
}
   10210:	90                   	nop
   10211:	c9                   	leave  
   10212:	c3                   	ret    
   10213:	90                   	nop

0000000000010214 <_ZNK3MPI3Win5FenceEi>:
}


inline void
MPI::Win::Fence(int assert) const
{
   10214:	f3 0f 1e fa          	endbr64 
   10218:	55                   	push   %rbp
   10219:	48 89 e5             	mov    %rsp,%rbp
   1021c:	48 83 ec 10          	sub    $0x10,%rsp
   10220:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10224:	89 75 f4             	mov    %esi,-0xc(%rbp)
  (void) MPI_Win_fence(assert, mpi_win);
   10227:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1022b:	48 8b 50 08          	mov    0x8(%rax),%rdx
   1022f:	8b 45 f4             	mov    -0xc(%rbp),%eax
   10232:	48 89 d6             	mov    %rdx,%rsi
   10235:	89 c7                	mov    %eax,%edi
   10237:	e8 24 a7 ff ff       	call   a960 <MPI_Win_fence@plt>
}
   1023c:	90                   	nop
   1023d:	c9                   	leave  
   1023e:	c3                   	ret    
   1023f:	90                   	nop

0000000000010240 <_ZNK3MPI3Win3GetEPKviRKNS_8DatatypeEiliS5_>:
MPI::Win::Get(const void *origin_addr, int origin_count,
		     const MPI::Datatype& origin_datatype,
		     int target_rank, MPI::Aint target_disp,
		     int target_count,
		     const MPI::Datatype& target_datatype) const
{
   10240:	f3 0f 1e fa          	endbr64 
   10244:	55                   	push   %rbp
   10245:	48 89 e5             	mov    %rsp,%rbp
   10248:	41 54                	push   %r12
   1024a:	53                   	push   %rbx
   1024b:	48 83 ec 30          	sub    $0x30,%rsp
   1024f:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   10253:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   10257:	89 55 dc             	mov    %edx,-0x24(%rbp)
   1025a:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
   1025e:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
   10262:	4c 89 4d c8          	mov    %r9,-0x38(%rbp)
  (void) MPI_Get(const_cast<void *>(origin_addr), origin_count, origin_datatype,
		     target_rank, target_disp,
		     target_count, target_datatype, mpi_win);
   10266:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1026a:	4c 8b 60 08          	mov    0x8(%rax),%r12
  (void) MPI_Get(const_cast<void *>(origin_addr), origin_count, origin_datatype,
   1026e:	48 8b 45 18          	mov    0x18(%rbp),%rax
   10272:	48 89 c7             	mov    %rax,%rdi
   10275:	e8 04 b6 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
   1027a:	48 89 c3             	mov    %rax,%rbx
   1027d:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   10281:	48 89 c7             	mov    %rax,%rdi
   10284:	e8 f5 b5 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
   10289:	48 89 c7             	mov    %rax,%rdi
   1028c:	48 8b 4d c8          	mov    -0x38(%rbp),%rcx
   10290:	8b 55 d8             	mov    -0x28(%rbp),%edx
   10293:	8b 75 dc             	mov    -0x24(%rbp),%esi
   10296:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   1029a:	41 54                	push   %r12
   1029c:	53                   	push   %rbx
   1029d:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
   102a1:	49 89 c8             	mov    %rcx,%r8
   102a4:	89 d1                	mov    %edx,%ecx
   102a6:	48 89 fa             	mov    %rdi,%rdx
   102a9:	48 89 c7             	mov    %rax,%rdi
   102ac:	e8 5f 9e ff ff       	call   a110 <MPI_Get@plt>
   102b1:	48 83 c4 10          	add    $0x10,%rsp

}
   102b5:	90                   	nop
   102b6:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
   102ba:	5b                   	pop    %rbx
   102bb:	41 5c                	pop    %r12
   102bd:	5d                   	pop    %rbp
   102be:	c3                   	ret    
   102bf:	90                   	nop

00000000000102c0 <_ZNK3MPI3Win9Get_groupEv>:


inline MPI::Group
MPI::Win::Get_group() const
{
   102c0:	f3 0f 1e fa          	endbr64 
   102c4:	55                   	push   %rbp
   102c5:	48 89 e5             	mov    %rsp,%rbp
   102c8:	48 83 ec 20          	sub    $0x20,%rsp
   102cc:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   102d0:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   102d4:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   102db:	00 00 
   102dd:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   102e1:	31 c0                	xor    %eax,%eax
  MPI_Group mpi_group;
  (void) MPI_Win_get_group(mpi_win, &mpi_group);
   102e3:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   102e7:	48 8b 40 08          	mov    0x8(%rax),%rax
   102eb:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
   102ef:	48 89 d6             	mov    %rdx,%rsi
   102f2:	48 89 c7             	mov    %rax,%rdi
   102f5:	e8 86 9e ff ff       	call   a180 <MPI_Win_get_group@plt>
  return mpi_group;
   102fa:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   102fe:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10302:	48 89 d6             	mov    %rdx,%rsi
   10305:	48 89 c7             	mov    %rax,%rdi
   10308:	e8 7f b7 ff ff       	call   ba8c <_ZN3MPI5GroupC1EP12ompi_group_t>
}
   1030d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10311:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   10318:	00 00 
   1031a:	74 05                	je     10321 <_ZNK3MPI3Win9Get_groupEv+0x61>
   1031c:	e8 bf a1 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   10321:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10325:	c9                   	leave  
   10326:	c3                   	ret    
   10327:	90                   	nop

0000000000010328 <_ZNK3MPI3Win4LockEiii>:


inline void
MPI::Win::Lock(int lock_type, int rank, int assert) const
{
   10328:	f3 0f 1e fa          	endbr64 
   1032c:	55                   	push   %rbp
   1032d:	48 89 e5             	mov    %rsp,%rbp
   10330:	48 83 ec 20          	sub    $0x20,%rsp
   10334:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10338:	89 75 f4             	mov    %esi,-0xc(%rbp)
   1033b:	89 55 f0             	mov    %edx,-0x10(%rbp)
   1033e:	89 4d ec             	mov    %ecx,-0x14(%rbp)
  (void) MPI_Win_lock(lock_type, rank, assert, mpi_win);
   10341:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10345:	48 8b 48 08          	mov    0x8(%rax),%rcx
   10349:	8b 55 ec             	mov    -0x14(%rbp),%edx
   1034c:	8b 75 f0             	mov    -0x10(%rbp),%esi
   1034f:	8b 45 f4             	mov    -0xc(%rbp),%eax
   10352:	89 c7                	mov    %eax,%edi
   10354:	e8 f7 a1 ff ff       	call   a550 <MPI_Win_lock@plt>
}
   10359:	90                   	nop
   1035a:	c9                   	leave  
   1035b:	c3                   	ret    

000000000001035c <_ZNK3MPI3Win4PostERKNS_5GroupEi>:


inline void
MPI::Win::Post(const MPI::Group& group, int assert) const
{
   1035c:	f3 0f 1e fa          	endbr64 
   10360:	55                   	push   %rbp
   10361:	48 89 e5             	mov    %rsp,%rbp
   10364:	53                   	push   %rbx
   10365:	48 83 ec 28          	sub    $0x28,%rsp
   10369:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   1036d:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   10371:	89 55 dc             	mov    %edx,-0x24(%rbp)
  (void) MPI_Win_post(group, assert, mpi_win);
   10374:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10378:	48 8b 58 08          	mov    0x8(%rax),%rbx
   1037c:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   10380:	48 89 c7             	mov    %rax,%rdi
   10383:	e8 80 b7 ff ff       	call   bb08 <_ZNK3MPI5GroupcvP12ompi_group_tEv>
   10388:	48 89 c1             	mov    %rax,%rcx
   1038b:	8b 45 dc             	mov    -0x24(%rbp),%eax
   1038e:	48 89 da             	mov    %rbx,%rdx
   10391:	89 c6                	mov    %eax,%esi
   10393:	48 89 cf             	mov    %rcx,%rdi
   10396:	e8 65 9c ff ff       	call   a000 <MPI_Win_post@plt>
}
   1039b:	90                   	nop
   1039c:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   103a0:	c9                   	leave  
   103a1:	c3                   	ret    

00000000000103a2 <_ZNK3MPI3Win3PutEPKviRKNS_8DatatypeEiliS5_>:
MPI::Win::Put(const void* origin_addr, int origin_count,
		     const MPI::Datatype& origin_datatype,
		     int target_rank, MPI::Aint target_disp,
		     int target_count,
		     const MPI::Datatype& target_datatype) const
{
   103a2:	f3 0f 1e fa          	endbr64 
   103a6:	55                   	push   %rbp
   103a7:	48 89 e5             	mov    %rsp,%rbp
   103aa:	41 54                	push   %r12
   103ac:	53                   	push   %rbx
   103ad:	48 83 ec 30          	sub    $0x30,%rsp
   103b1:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   103b5:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   103b9:	89 55 dc             	mov    %edx,-0x24(%rbp)
   103bc:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
   103c0:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
   103c4:	4c 89 4d c8          	mov    %r9,-0x38(%rbp)
  (void) MPI_Put(const_cast<void *>(origin_addr), origin_count, origin_datatype,
		 target_rank, target_disp, target_count,
		 target_datatype, mpi_win);
   103c8:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   103cc:	4c 8b 60 08          	mov    0x8(%rax),%r12
  (void) MPI_Put(const_cast<void *>(origin_addr), origin_count, origin_datatype,
   103d0:	48 8b 45 18          	mov    0x18(%rbp),%rax
   103d4:	48 89 c7             	mov    %rax,%rdi
   103d7:	e8 a2 b4 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
   103dc:	48 89 c3             	mov    %rax,%rbx
   103df:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   103e3:	48 89 c7             	mov    %rax,%rdi
   103e6:	e8 93 b4 ff ff       	call   b87e <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
   103eb:	48 89 c7             	mov    %rax,%rdi
   103ee:	48 8b 4d c8          	mov    -0x38(%rbp),%rcx
   103f2:	8b 55 d8             	mov    -0x28(%rbp),%edx
   103f5:	8b 75 dc             	mov    -0x24(%rbp),%esi
   103f8:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   103fc:	41 54                	push   %r12
   103fe:	53                   	push   %rbx
   103ff:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
   10403:	49 89 c8             	mov    %rcx,%r8
   10406:	89 d1                	mov    %edx,%ecx
   10408:	48 89 fa             	mov    %rdi,%rdx
   1040b:	48 89 c7             	mov    %rax,%rdi
   1040e:	e8 ad a1 ff ff       	call   a5c0 <MPI_Put@plt>
   10413:	48 83 c4 10          	add    $0x10,%rsp

}
   10417:	90                   	nop
   10418:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
   1041c:	5b                   	pop    %rbx
   1041d:	41 5c                	pop    %r12
   1041f:	5d                   	pop    %rbp
   10420:	c3                   	ret    
   10421:	90                   	nop

0000000000010422 <_ZNK3MPI3Win5StartERKNS_5GroupEi>:


inline void
MPI::Win::Start(const MPI::Group& group, int assert) const
{
   10422:	f3 0f 1e fa          	endbr64 
   10426:	55                   	push   %rbp
   10427:	48 89 e5             	mov    %rsp,%rbp
   1042a:	53                   	push   %rbx
   1042b:	48 83 ec 28          	sub    $0x28,%rsp
   1042f:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   10433:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   10437:	89 55 dc             	mov    %edx,-0x24(%rbp)
  (void) MPI_Win_start(group, assert, mpi_win);
   1043a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1043e:	48 8b 58 08          	mov    0x8(%rax),%rbx
   10442:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   10446:	48 89 c7             	mov    %rax,%rdi
   10449:	e8 ba b6 ff ff       	call   bb08 <_ZNK3MPI5GroupcvP12ompi_group_tEv>
   1044e:	48 89 c1             	mov    %rax,%rcx
   10451:	8b 45 dc             	mov    -0x24(%rbp),%eax
   10454:	48 89 da             	mov    %rbx,%rdx
   10457:	89 c6                	mov    %eax,%esi
   10459:	48 89 cf             	mov    %rcx,%rdi
   1045c:	e8 8f 9d ff ff       	call   a1f0 <MPI_Win_start@plt>
}
   10461:	90                   	nop
   10462:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   10466:	c9                   	leave  
   10467:	c3                   	ret    

0000000000010468 <_ZNK3MPI3Win4TestEv>:


inline bool
MPI::Win::Test() const
{
   10468:	f3 0f 1e fa          	endbr64 
   1046c:	55                   	push   %rbp
   1046d:	48 89 e5             	mov    %rsp,%rbp
   10470:	48 83 ec 20          	sub    $0x20,%rsp
   10474:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   10478:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   1047f:	00 00 
   10481:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   10485:	31 c0                	xor    %eax,%eax
  int flag;
  MPI_Win_test(mpi_win, &flag);
   10487:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1048b:	48 8b 40 08          	mov    0x8(%rax),%rax
   1048f:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
   10493:	48 89 d6             	mov    %rdx,%rsi
   10496:	48 89 c7             	mov    %rax,%rdi
   10499:	e8 62 a4 ff ff       	call   a900 <MPI_Win_test@plt>
  return OPAL_INT_TO_BOOL(flag);
   1049e:	8b 45 f4             	mov    -0xc(%rbp),%eax
   104a1:	85 c0                	test   %eax,%eax
   104a3:	0f 95 c0             	setne  %al

}
   104a6:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   104aa:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
   104b1:	00 00 
   104b3:	74 05                	je     104ba <_ZNK3MPI3Win4TestEv+0x52>
   104b5:	e8 26 a0 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   104ba:	c9                   	leave  
   104bb:	c3                   	ret    

00000000000104bc <_ZNK3MPI3Win6UnlockEi>:


inline void
MPI::Win::Unlock(int rank) const
{
   104bc:	f3 0f 1e fa          	endbr64 
   104c0:	55                   	push   %rbp
   104c1:	48 89 e5             	mov    %rsp,%rbp
   104c4:	48 83 ec 10          	sub    $0x10,%rsp
   104c8:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   104cc:	89 75 f4             	mov    %esi,-0xc(%rbp)
  (void) MPI_Win_unlock(rank, mpi_win);
   104cf:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   104d3:	48 8b 50 08          	mov    0x8(%rax),%rdx
   104d7:	8b 45 f4             	mov    -0xc(%rbp),%eax
   104da:	48 89 d6             	mov    %rdx,%rsi
   104dd:	89 c7                	mov    %eax,%edi
   104df:	e8 0c 99 ff ff       	call   9df0 <MPI_Win_unlock@plt>
}
   104e4:	90                   	nop
   104e5:	c9                   	leave  
   104e6:	c3                   	ret    
   104e7:	90                   	nop

00000000000104e8 <_ZNK3MPI3Win4WaitEv>:


inline void
MPI::Win::Wait() const
{
   104e8:	f3 0f 1e fa          	endbr64 
   104ec:	55                   	push   %rbp
   104ed:	48 89 e5             	mov    %rsp,%rbp
   104f0:	48 83 ec 10          	sub    $0x10,%rsp
   104f4:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
  (void) MPI_Win_wait(mpi_win);
   104f8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   104fc:	48 8b 40 08          	mov    0x8(%rax),%rax
   10500:	48 89 c7             	mov    %rax,%rdi
   10503:	e8 b8 9d ff ff       	call   a2c0 <MPI_Win_wait@plt>
}
   10508:	90                   	nop
   10509:	c9                   	leave  
   1050a:	c3                   	ret    
   1050b:	90                   	nop

000000000001050c <_ZNK3MPI3Win15Call_errhandlerEi>:
// External Interfaces
//

inline void
MPI::Win::Call_errhandler(int errorcode) const
{
   1050c:	f3 0f 1e fa          	endbr64 
   10510:	55                   	push   %rbp
   10511:	48 89 e5             	mov    %rsp,%rbp
   10514:	48 83 ec 10          	sub    $0x10,%rsp
   10518:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   1051c:	89 75 f4             	mov    %esi,-0xc(%rbp)
  (void) MPI_Win_call_errhandler(mpi_win, errorcode);
   1051f:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10523:	48 8b 40 08          	mov    0x8(%rax),%rax
   10527:	8b 55 f4             	mov    -0xc(%rbp),%edx
   1052a:	89 d6                	mov    %edx,%esi
   1052c:	48 89 c7             	mov    %rax,%rdi
   1052f:	e8 fc a3 ff ff       	call   a930 <MPI_Win_call_errhandler@plt>
}
   10534:	90                   	nop
   10535:	c9                   	leave  
   10536:	c3                   	ret    
   10537:	90                   	nop

0000000000010538 <_ZN3MPI3Win11Delete_attrEi>:
    return (MPI_SUCCESS == ret) ? keyval : ret;
}

inline void
MPI::Win::Delete_attr(int win_keyval)
{
   10538:	f3 0f 1e fa          	endbr64 
   1053c:	55                   	push   %rbp
   1053d:	48 89 e5             	mov    %rsp,%rbp
   10540:	48 83 ec 10          	sub    $0x10,%rsp
   10544:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10548:	89 75 f4             	mov    %esi,-0xc(%rbp)
  (void) MPI_Win_delete_attr(mpi_win, win_keyval);
   1054b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1054f:	48 8b 40 08          	mov    0x8(%rax),%rax
   10553:	8b 55 f4             	mov    -0xc(%rbp),%edx
   10556:	89 d6                	mov    %edx,%esi
   10558:	48 89 c7             	mov    %rax,%rdi
   1055b:	e8 a0 a4 ff ff       	call   aa00 <MPI_Win_delete_attr@plt>
}
   10560:	90                   	nop
   10561:	c9                   	leave  
   10562:	c3                   	ret    
   10563:	90                   	nop

0000000000010564 <_ZNK3MPI3Win8Get_nameEPcRi>:
}


inline void
MPI::Win::Get_name(char* win_name, int& resultlen) const
{
   10564:	f3 0f 1e fa          	endbr64 
   10568:	55                   	push   %rbp
   10569:	48 89 e5             	mov    %rsp,%rbp
   1056c:	48 83 ec 20          	sub    $0x20,%rsp
   10570:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10574:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   10578:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
  (void) MPI_Win_get_name(mpi_win, win_name, &resultlen);
   1057c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10580:	48 8b 40 08          	mov    0x8(%rax),%rax
   10584:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   10588:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
   1058c:	48 89 ce             	mov    %rcx,%rsi
   1058f:	48 89 c7             	mov    %rax,%rdi
   10592:	e8 79 99 ff ff       	call   9f10 <MPI_Win_get_name@plt>
}
   10597:	90                   	nop
   10598:	c9                   	leave  
   10599:	c3                   	ret    

000000000001059a <_ZN3MPI3Win8Set_attrEiPKv>:


inline void
MPI::Win::Set_attr(int win_keyval, const void* attribute_val)
{
   1059a:	f3 0f 1e fa          	endbr64 
   1059e:	55                   	push   %rbp
   1059f:	48 89 e5             	mov    %rsp,%rbp
   105a2:	48 83 ec 20          	sub    $0x20,%rsp
   105a6:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   105aa:	89 75 f4             	mov    %esi,-0xc(%rbp)
   105ad:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
  (void) MPI_Win_set_attr(mpi_win, win_keyval, const_cast<void *>(attribute_val));
   105b1:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   105b5:	48 8b 40 08          	mov    0x8(%rax),%rax
   105b9:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   105bd:	8b 4d f4             	mov    -0xc(%rbp),%ecx
   105c0:	89 ce                	mov    %ecx,%esi
   105c2:	48 89 c7             	mov    %rax,%rdi
   105c5:	e8 26 a4 ff ff       	call   a9f0 <MPI_Win_set_attr@plt>
}
   105ca:	90                   	nop
   105cb:	c9                   	leave  
   105cc:	c3                   	ret    
   105cd:	90                   	nop

00000000000105ce <_ZN3MPI3Win8Set_nameEPKc>:


inline void
MPI::Win::Set_name(const char* win_name)
{
   105ce:	f3 0f 1e fa          	endbr64 
   105d2:	55                   	push   %rbp
   105d3:	48 89 e5             	mov    %rsp,%rbp
   105d6:	48 83 ec 10          	sub    $0x10,%rsp
   105da:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   105de:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
  (void) MPI_Win_set_name(mpi_win, const_cast<char *>(win_name));
   105e2:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   105e6:	48 8b 40 08          	mov    0x8(%rax),%rax
   105ea:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   105ee:	48 89 d6             	mov    %rdx,%rsi
   105f1:	48 89 c7             	mov    %rax,%rdi
   105f4:	e8 77 9d ff ff       	call   a370 <MPI_Win_set_name@plt>
}
   105f9:	90                   	nop
   105fa:	c9                   	leave  
   105fb:	c3                   	ret    

00000000000105fc <_ZNSt8__detail18__to_chars_10_implIjEEvPcjT_>:
  // Write an unsigned integer value to the range [first,first+len).
  // The caller is required to provide a buffer of exactly the right size
  // (which can be determined by the __to_chars_len function).
  template<typename _Tp>
    void
    __to_chars_10_impl(char* __first, unsigned __len, _Tp __val) noexcept
   105fc:	f3 0f 1e fa          	endbr64 
   10600:	55                   	push   %rbp
   10601:	48 89 e5             	mov    %rsp,%rbp
   10604:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   10608:	89 75 e4             	mov    %esi,-0x1c(%rbp)
   1060b:	89 55 e0             	mov    %edx,-0x20(%rbp)
	"0001020304050607080910111213141516171819"
	"2021222324252627282930313233343536373839"
	"4041424344454647484950515253545556575859"
	"6061626364656667686970717273747576777879"
	"8081828384858687888990919293949596979899";
      unsigned __pos = __len - 1;
   1060e:	8b 45 e4             	mov    -0x1c(%rbp),%eax
   10611:	83 e8 01             	sub    $0x1,%eax
   10614:	89 45 f4             	mov    %eax,-0xc(%rbp)
      while (__val >= 100)
   10617:	eb 77                	jmp    10690 <_ZNSt8__detail18__to_chars_10_implIjEEvPcjT_+0x94>
	{
	  auto const __num = (__val % 100) * 2;
   10619:	8b 55 e0             	mov    -0x20(%rbp),%edx
   1061c:	89 d0                	mov    %edx,%eax
   1061e:	48 69 c0 1f 85 eb 51 	imul   $0x51eb851f,%rax,%rax
   10625:	48 c1 e8 20          	shr    $0x20,%rax
   10629:	c1 e8 05             	shr    $0x5,%eax
   1062c:	6b c8 64             	imul   $0x64,%eax,%ecx
   1062f:	89 d0                	mov    %edx,%eax
   10631:	29 c8                	sub    %ecx,%eax
   10633:	01 c0                	add    %eax,%eax
   10635:	89 45 fc             	mov    %eax,-0x4(%rbp)
	  __val /= 100;
   10638:	8b 45 e0             	mov    -0x20(%rbp),%eax
   1063b:	89 c0                	mov    %eax,%eax
   1063d:	48 69 c0 1f 85 eb 51 	imul   $0x51eb851f,%rax,%rax
   10644:	48 c1 e8 20          	shr    $0x20,%rax
   10648:	c1 e8 05             	shr    $0x5,%eax
   1064b:	89 45 e0             	mov    %eax,-0x20(%rbp)
	  __first[__pos] = __digits[__num + 1];
   1064e:	8b 45 fc             	mov    -0x4(%rbp),%eax
   10651:	8d 48 01             	lea    0x1(%rax),%ecx
   10654:	8b 55 f4             	mov    -0xc(%rbp),%edx
   10657:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1065b:	48 01 c2             	add    %rax,%rdx
   1065e:	89 c9                	mov    %ecx,%ecx
   10660:	48 8d 05 d9 5e 00 00 	lea    0x5ed9(%rip),%rax        # 16540 <_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits>
   10667:	0f b6 04 01          	movzbl (%rcx,%rax,1),%eax
   1066b:	88 02                	mov    %al,(%rdx)
	  __first[__pos - 1] = __digits[__num];
   1066d:	8b 45 f4             	mov    -0xc(%rbp),%eax
   10670:	83 e8 01             	sub    $0x1,%eax
   10673:	89 c2                	mov    %eax,%edx
   10675:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10679:	48 01 c2             	add    %rax,%rdx
   1067c:	8b 45 fc             	mov    -0x4(%rbp),%eax
   1067f:	48 8d 0d ba 5e 00 00 	lea    0x5eba(%rip),%rcx        # 16540 <_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits>
   10686:	0f b6 04 08          	movzbl (%rax,%rcx,1),%eax
   1068a:	88 02                	mov    %al,(%rdx)
	  __pos -= 2;
   1068c:	83 6d f4 02          	subl   $0x2,-0xc(%rbp)
      while (__val >= 100)
   10690:	83 7d e0 63          	cmpl   $0x63,-0x20(%rbp)
   10694:	77 83                	ja     10619 <_ZNSt8__detail18__to_chars_10_implIjEEvPcjT_+0x1d>
	}
      if (__val >= 10)
   10696:	83 7d e0 09          	cmpl   $0x9,-0x20(%rbp)
   1069a:	76 3b                	jbe    106d7 <_ZNSt8__detail18__to_chars_10_implIjEEvPcjT_+0xdb>
	{
	  auto const __num = __val * 2;
   1069c:	8b 45 e0             	mov    -0x20(%rbp),%eax
   1069f:	01 c0                	add    %eax,%eax
   106a1:	89 45 f8             	mov    %eax,-0x8(%rbp)
	  __first[1] = __digits[__num + 1];
   106a4:	8b 45 f8             	mov    -0x8(%rbp),%eax
   106a7:	8d 48 01             	lea    0x1(%rax),%ecx
   106aa:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   106ae:	48 8d 50 01          	lea    0x1(%rax),%rdx
   106b2:	89 c9                	mov    %ecx,%ecx
   106b4:	48 8d 05 85 5e 00 00 	lea    0x5e85(%rip),%rax        # 16540 <_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits>
   106bb:	0f b6 04 01          	movzbl (%rcx,%rax,1),%eax
   106bf:	88 02                	mov    %al,(%rdx)
	  __first[0] = __digits[__num];
   106c1:	8b 45 f8             	mov    -0x8(%rbp),%eax
   106c4:	48 8d 15 75 5e 00 00 	lea    0x5e75(%rip),%rdx        # 16540 <_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits>
   106cb:	0f b6 14 10          	movzbl (%rax,%rdx,1),%edx
   106cf:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   106d3:	88 10                	mov    %dl,(%rax)
	}
      else
	__first[0] = '0' + __val;
    }
   106d5:	eb 0e                	jmp    106e5 <_ZNSt8__detail18__to_chars_10_implIjEEvPcjT_+0xe9>
	__first[0] = '0' + __val;
   106d7:	8b 45 e0             	mov    -0x20(%rbp),%eax
   106da:	83 c0 30             	add    $0x30,%eax
   106dd:	89 c2                	mov    %eax,%edx
   106df:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   106e3:	88 10                	mov    %dl,(%rax)
    }
   106e5:	90                   	nop
   106e6:	5d                   	pop    %rbp
   106e7:	c3                   	ret    

00000000000106e8 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_>:
    operator+(basic_string<_CharT, _Traits, _Alloc>&& __lhs,
   106e8:	f3 0f 1e fa          	endbr64 
   106ec:	55                   	push   %rbp
   106ed:	48 89 e5             	mov    %rsp,%rbp
   106f0:	53                   	push   %rbx
   106f1:	48 83 ec 48          	sub    $0x48,%rsp
   106f5:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
   106f9:	48 89 75 c0          	mov    %rsi,-0x40(%rbp)
   106fd:	48 89 55 b8          	mov    %rdx,-0x48(%rbp)
   10701:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   10708:	00 00 
   1070a:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   1070e:	31 c0                	xor    %eax,%eax
      bool __use_rhs = false;
   10710:	c6 45 df 00          	movb   $0x0,-0x21(%rbp)
      if _GLIBCXX17_CONSTEXPR (typename _Alloc_traits::is_always_equal{})
   10714:	48 8d 45 de          	lea    -0x22(%rbp),%rax
   10718:	48 89 c7             	mov    %rax,%rdi
   1071b:	e8 a8 b0 ff ff       	call   b7c8 <_ZNKSt17integral_constantIbLb1EEcvbEv>
   10720:	84 c0                	test   %al,%al
   10722:	74 06                	je     1072a <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x42>
	__use_rhs = true;
   10724:	c6 45 df 01          	movb   $0x1,-0x21(%rbp)
   10728:	eb 5b                	jmp    10785 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x9d>
      else if (__lhs.get_allocator() == __rhs.get_allocator())
   1072a:	48 8d 45 de          	lea    -0x22(%rbp),%rax
   1072e:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
   10732:	48 89 d6             	mov    %rdx,%rsi
   10735:	48 89 c7             	mov    %rax,%rdi
   10738:	e8 73 98 ff ff       	call   9fb0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE13get_allocatorEv@plt>
   1073d:	48 8d 45 dd          	lea    -0x23(%rbp),%rax
   10741:	48 8b 55 c0          	mov    -0x40(%rbp),%rdx
   10745:	48 89 d6             	mov    %rdx,%rsi
   10748:	48 89 c7             	mov    %rax,%rdi
   1074b:	e8 60 98 ff ff       	call   9fb0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE13get_allocatorEv@plt>
   10750:	48 8d 55 de          	lea    -0x22(%rbp),%rdx
   10754:	48 8d 45 dd          	lea    -0x23(%rbp),%rax
   10758:	48 89 d6             	mov    %rdx,%rsi
   1075b:	48 89 c7             	mov    %rax,%rdi
   1075e:	e8 4b 01 00 00       	call   108ae <_ZSteqRKSaIcES1_>
   10763:	89 c3                	mov    %eax,%ebx
   10765:	48 8d 45 dd          	lea    -0x23(%rbp),%rax
   10769:	48 89 c7             	mov    %rax,%rdi
   1076c:	e8 ff 9c ff ff       	call   a470 <_ZNSaIcED1Ev@plt>
   10771:	48 8d 45 de          	lea    -0x22(%rbp),%rax
   10775:	48 89 c7             	mov    %rax,%rdi
   10778:	e8 f3 9c ff ff       	call   a470 <_ZNSaIcED1Ev@plt>
   1077d:	84 db                	test   %bl,%bl
   1077f:	74 04                	je     10785 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x9d>
	__use_rhs = true;
   10781:	c6 45 df 01          	movb   $0x1,-0x21(%rbp)
      if (__use_rhs)
   10785:	80 7d df 00          	cmpb   $0x0,-0x21(%rbp)
   10789:	0f 84 87 00 00 00    	je     10816 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x12e>
	  const auto __size = __lhs.size() + __rhs.size();
   1078f:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   10793:	48 89 c7             	mov    %rax,%rdi
   10796:	e8 15 9a ff ff       	call   a1b0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE4sizeEv@plt>
   1079b:	48 89 c3             	mov    %rax,%rbx
   1079e:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   107a2:	48 89 c7             	mov    %rax,%rdi
   107a5:	e8 06 9a ff ff       	call   a1b0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE4sizeEv@plt>
   107aa:	48 01 d8             	add    %rbx,%rax
   107ad:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
	  if (__size > __lhs.capacity() && __size <= __rhs.capacity())
   107b1:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   107b5:	48 89 c7             	mov    %rax,%rdi
   107b8:	e8 13 96 ff ff       	call   9dd0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE8capacityEv@plt>
   107bd:	48 39 45 e0          	cmp    %rax,-0x20(%rbp)
   107c1:	76 19                	jbe    107dc <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0xf4>
   107c3:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   107c7:	48 89 c7             	mov    %rax,%rdi
   107ca:	e8 01 96 ff ff       	call   9dd0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE8capacityEv@plt>
   107cf:	48 39 45 e0          	cmp    %rax,-0x20(%rbp)
   107d3:	77 07                	ja     107dc <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0xf4>
   107d5:	b8 01 00 00 00       	mov    $0x1,%eax
   107da:	eb 05                	jmp    107e1 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0xf9>
   107dc:	b8 00 00 00 00       	mov    $0x0,%eax
   107e1:	84 c0                	test   %al,%al
   107e3:	74 31                	je     10816 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x12e>
	    return std::move(__rhs.insert(0, __lhs));
   107e5:	48 8b 55 c0          	mov    -0x40(%rbp),%rdx
   107e9:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   107ed:	be 00 00 00 00       	mov    $0x0,%esi
   107f2:	48 89 c7             	mov    %rax,%rdi
   107f5:	e8 56 a0 ff ff       	call   a850 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmRKS4_@plt>
   107fa:	48 89 c7             	mov    %rax,%rdi
   107fd:	e8 c3 00 00 00       	call   108c5 <_ZSt4moveIRNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEONSt16remove_referenceIT_E4typeEOS8_>
   10802:	48 89 c2             	mov    %rax,%rdx
   10805:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   10809:	48 89 d6             	mov    %rdx,%rsi
   1080c:	48 89 c7             	mov    %rax,%rdi
   1080f:	e8 dc 98 ff ff       	call   a0f0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@plt>
   10814:	eb 2d                	jmp    10843 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x15b>
      return std::move(__lhs.append(__rhs));
   10816:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
   1081a:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   1081e:	48 89 d6             	mov    %rdx,%rsi
   10821:	48 89 c7             	mov    %rax,%rdi
   10824:	e8 d7 a2 ff ff       	call   ab00 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendERKS4_@plt>
   10829:	48 89 c7             	mov    %rax,%rdi
   1082c:	e8 94 00 00 00       	call   108c5 <_ZSt4moveIRNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEONSt16remove_referenceIT_E4typeEOS8_>
   10831:	48 89 c2             	mov    %rax,%rdx
   10834:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   10838:	48 89 d6             	mov    %rdx,%rsi
   1083b:	48 89 c7             	mov    %rax,%rdi
   1083e:	e8 ad 98 ff ff       	call   a0f0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@plt>
    }
   10843:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10847:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   1084e:	00 00 
   10850:	74 05                	je     10857 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x16f>
   10852:	e8 89 9c ff ff       	call   a4e0 <__stack_chk_fail@plt>
   10857:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   1085b:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   1085f:	c9                   	leave  
   10860:	c3                   	ret    

0000000000010861 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEPKS5_OS8_>:
    operator+(const _CharT* __lhs,
   10861:	f3 0f 1e fa          	endbr64 
   10865:	55                   	push   %rbp
   10866:	48 89 e5             	mov    %rsp,%rbp
   10869:	48 83 ec 20          	sub    $0x20,%rsp
   1086d:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10871:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   10875:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    { return std::move(__rhs.insert(0, __lhs)); }
   10879:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   1087d:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10881:	be 00 00 00 00       	mov    $0x0,%esi
   10886:	48 89 c7             	mov    %rax,%rdi
   10889:	e8 72 96 ff ff       	call   9f00 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmPKc@plt>
   1088e:	48 89 c7             	mov    %rax,%rdi
   10891:	e8 2f 00 00 00       	call   108c5 <_ZSt4moveIRNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEONSt16remove_referenceIT_E4typeEOS8_>
   10896:	48 89 c2             	mov    %rax,%rdx
   10899:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1089d:	48 89 d6             	mov    %rdx,%rsi
   108a0:	48 89 c7             	mov    %rax,%rdi
   108a3:	e8 48 98 ff ff       	call   a0f0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@plt>
   108a8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   108ac:	c9                   	leave  
   108ad:	c3                   	ret    

00000000000108ae <_ZSteqRKSaIcES1_>:
	__allocator_base<_Tp>::deallocate(__p, __n);
      }
#endif // C++20

      friend _GLIBCXX20_CONSTEXPR bool
      operator==(const allocator&, const allocator&) _GLIBCXX_NOTHROW
   108ae:	f3 0f 1e fa          	endbr64 
   108b2:	55                   	push   %rbp
   108b3:	48 89 e5             	mov    %rsp,%rbp
   108b6:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   108ba:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
      { return true; }
   108be:	b8 01 00 00 00       	mov    $0x1,%eax
   108c3:	5d                   	pop    %rbp
   108c4:	c3                   	ret    

00000000000108c5 <_ZSt4moveIRNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEONSt16remove_referenceIT_E4typeEOS8_>:
   *  @return The parameter cast to an rvalue-reference to allow moving it.
  */
  template<typename _Tp>
    _GLIBCXX_NODISCARD
    constexpr typename std::remove_reference<_Tp>::type&&
    move(_Tp&& __t) noexcept
   108c5:	f3 0f 1e fa          	endbr64 
   108c9:	55                   	push   %rbp
   108ca:	48 89 e5             	mov    %rsp,%rbp
   108cd:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    { return static_cast<typename std::remove_reference<_Tp>::type&&>(__t); }
   108d1:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   108d5:	5d                   	pop    %rbp
   108d6:	c3                   	ret    
   108d7:	90                   	nop

00000000000108d8 <_ZN3MPI9IntercommD1Ev>:
class Intercomm : public Comm {
   108d8:	f3 0f 1e fa          	endbr64 
   108dc:	55                   	push   %rbp
   108dd:	48 89 e5             	mov    %rsp,%rbp
   108e0:	48 83 ec 10          	sub    $0x10,%rsp
   108e4:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   108e8:	48 8d 15 a1 ae 00 00 	lea    0xaea1(%rip),%rdx        # 1b790 <_ZTVN3MPI9IntercommE+0x10>
   108ef:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   108f3:	48 89 10             	mov    %rdx,(%rax)
   108f6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   108fa:	48 89 c7             	mov    %rax,%rdi
   108fd:	e8 c6 b3 ff ff       	call   bcc8 <_ZN3MPI4CommD1Ev>
   10902:	90                   	nop
   10903:	c9                   	leave  
   10904:	c3                   	ret    
   10905:	90                   	nop

0000000000010906 <_ZN3MPI9IntercommD0Ev>:
   10906:	f3 0f 1e fa          	endbr64 
   1090a:	55                   	push   %rbp
   1090b:	48 89 e5             	mov    %rsp,%rbp
   1090e:	48 83 ec 10          	sub    $0x10,%rsp
   10912:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10916:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1091a:	48 89 c7             	mov    %rax,%rdi
   1091d:	e8 b6 ff ff ff       	call   108d8 <_ZN3MPI9IntercommD1Ev>
   10922:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10926:	be 10 00 00 00       	mov    $0x10,%esi
   1092b:	48 89 c7             	mov    %rax,%rdi
   1092e:	e8 cd 9a ff ff       	call   a400 <_ZdlPvm@plt>
   10933:	c9                   	leave  
   10934:	c3                   	ret    
   10935:	90                   	nop

0000000000010936 <_ZN3MPI9GraphcommD1Ev>:

//===================================================================
//                    Class Graphcomm
//===================================================================

class Graphcomm : public Intracomm {
   10936:	f3 0f 1e fa          	endbr64 
   1093a:	55                   	push   %rbp
   1093b:	48 89 e5             	mov    %rsp,%rbp
   1093e:	48 83 ec 10          	sub    $0x10,%rsp
   10942:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10946:	48 8d 15 3b b0 00 00 	lea    0xb03b(%rip),%rdx        # 1b988 <_ZTVN3MPI9GraphcommE+0x10>
   1094d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10951:	48 89 10             	mov    %rdx,(%rax)
   10954:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10958:	48 89 c7             	mov    %rax,%rdi
   1095b:	e8 f4 b3 ff ff       	call   bd54 <_ZN3MPI9IntracommD1Ev>
   10960:	90                   	nop
   10961:	c9                   	leave  
   10962:	c3                   	ret    
   10963:	90                   	nop

0000000000010964 <_ZN3MPI9GraphcommD0Ev>:
   10964:	f3 0f 1e fa          	endbr64 
   10968:	55                   	push   %rbp
   10969:	48 89 e5             	mov    %rsp,%rbp
   1096c:	48 83 ec 10          	sub    $0x10,%rsp
   10970:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10974:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10978:	48 89 c7             	mov    %rax,%rdi
   1097b:	e8 b6 ff ff ff       	call   10936 <_ZN3MPI9GraphcommD1Ev>
   10980:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10984:	be 10 00 00 00       	mov    $0x10,%esi
   10989:	48 89 c7             	mov    %rax,%rdi
   1098c:	e8 6f 9a ff ff       	call   a400 <_ZdlPvm@plt>
   10991:	c9                   	leave  
   10992:	c3                   	ret    
   10993:	90                   	nop

0000000000010994 <_ZN3MPI8CartcommD1Ev>:
class Cartcomm : public Intracomm {
   10994:	f3 0f 1e fa          	endbr64 
   10998:	55                   	push   %rbp
   10999:	48 89 e5             	mov    %rsp,%rbp
   1099c:	48 83 ec 10          	sub    $0x10,%rsp
   109a0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   109a4:	48 8d 15 3d b2 00 00 	lea    0xb23d(%rip),%rdx        # 1bbe8 <_ZTVN3MPI8CartcommE+0x10>
   109ab:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   109af:	48 89 10             	mov    %rdx,(%rax)
   109b2:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   109b6:	48 89 c7             	mov    %rax,%rdi
   109b9:	e8 96 b3 ff ff       	call   bd54 <_ZN3MPI9IntracommD1Ev>
   109be:	90                   	nop
   109bf:	c9                   	leave  
   109c0:	c3                   	ret    
   109c1:	90                   	nop

00000000000109c2 <_ZN3MPI8CartcommD0Ev>:
   109c2:	f3 0f 1e fa          	endbr64 
   109c6:	55                   	push   %rbp
   109c7:	48 89 e5             	mov    %rsp,%rbp
   109ca:	48 83 ec 10          	sub    $0x10,%rsp
   109ce:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   109d2:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   109d6:	48 89 c7             	mov    %rax,%rdi
   109d9:	e8 b6 ff ff ff       	call   10994 <_ZN3MPI8CartcommD1Ev>
   109de:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   109e2:	be 10 00 00 00       	mov    $0x10,%esi
   109e7:	48 89 c7             	mov    %rax,%rdi
   109ea:	e8 11 9a ff ff       	call   a400 <_ZdlPvm@plt>
   109ef:	c9                   	leave  
   109f0:	c3                   	ret    

00000000000109f1 <_Z14updateVelocityP10particle_td>:
#include <stdlib.h>
#include "string.h"
#include "stdio.h"


void updateVelocity(particle_t* particle, double timeStep) {
   109f1:	f3 0f 1e fa          	endbr64 
   109f5:	55                   	push   %rbp
   109f6:	48 89 e5             	mov    %rsp,%rbp
   109f9:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   109fd:	f2 0f 11 45 d0       	movsd  %xmm0,-0x30(%rbp)
    if(!particle->enabled)
   10a02:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10a06:	0f b6 40 38          	movzbl 0x38(%rax),%eax
   10a0a:	83 f0 01             	xor    $0x1,%eax
   10a0d:	84 c0                	test   %al,%al
   10a0f:	75 64                	jne    10a75 <_Z14updateVelocityP10particle_td+0x84>
        return;
    double* accelerationDim = particle->acceleration;
   10a11:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10a15:	48 8b 40 18          	mov    0x18(%rax),%rax
   10a19:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
    double* velocityDim = particle->velocity;
   10a1d:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10a21:	48 8b 40 10          	mov    0x10(%rax),%rax
   10a25:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    #pragma unroll
    for(int n=0; n<particle->ndim; n++) 
   10a29:	c7 45 ec 00 00 00 00 	movl   $0x0,-0x14(%rbp)
   10a30:	eb 35                	jmp    10a67 <_Z14updateVelocityP10particle_td+0x76>
        *(velocityDim++) += *(accelerationDim++) * timeStep;
   10a32:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   10a36:	48 8d 50 08          	lea    0x8(%rax),%rdx
   10a3a:	48 89 55 f0          	mov    %rdx,-0x10(%rbp)
   10a3e:	f2 0f 10 00          	movsd  (%rax),%xmm0
   10a42:	66 0f 28 c8          	movapd %xmm0,%xmm1
   10a46:	f2 0f 59 4d d0       	mulsd  -0x30(%rbp),%xmm1
   10a4b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10a4f:	48 8d 50 08          	lea    0x8(%rax),%rdx
   10a53:	48 89 55 f8          	mov    %rdx,-0x8(%rbp)
   10a57:	f2 0f 10 00          	movsd  (%rax),%xmm0
   10a5b:	f2 0f 58 c1          	addsd  %xmm1,%xmm0
   10a5f:	f2 0f 11 00          	movsd  %xmm0,(%rax)
    for(int n=0; n<particle->ndim; n++) 
   10a63:	83 45 ec 01          	addl   $0x1,-0x14(%rbp)
   10a67:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10a6b:	8b 40 30             	mov    0x30(%rax),%eax
   10a6e:	39 45 ec             	cmp    %eax,-0x14(%rbp)
   10a71:	7c bf                	jl     10a32 <_Z14updateVelocityP10particle_td+0x41>
   10a73:	eb 01                	jmp    10a76 <_Z14updateVelocityP10particle_td+0x85>
        return;
   10a75:	90                   	nop
}
   10a76:	5d                   	pop    %rbp
   10a77:	c3                   	ret    

0000000000010a78 <_Z14updatePositionP10particle_td>:

void updatePosition(particle_t* particle, double timeStep) {
   10a78:	f3 0f 1e fa          	endbr64 
   10a7c:	55                   	push   %rbp
   10a7d:	48 89 e5             	mov    %rsp,%rbp
   10a80:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   10a84:	f2 0f 11 45 d0       	movsd  %xmm0,-0x30(%rbp)
    if(!particle->enabled)
   10a89:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10a8d:	0f b6 40 38          	movzbl 0x38(%rax),%eax
   10a91:	83 f0 01             	xor    $0x1,%eax
   10a94:	84 c0                	test   %al,%al
   10a96:	75 64                	jne    10afc <_Z14updatePositionP10particle_td+0x84>
        return;
    double* positionDim = particle->position;
   10a98:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10a9c:	48 8b 40 08          	mov    0x8(%rax),%rax
   10aa0:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
    double* velocityDim = particle->velocity;
   10aa4:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10aa8:	48 8b 40 10          	mov    0x10(%rax),%rax
   10aac:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    #pragma unroll
    for(int n=0; n<particle->ndim; n++) 
   10ab0:	c7 45 ec 00 00 00 00 	movl   $0x0,-0x14(%rbp)
   10ab7:	eb 35                	jmp    10aee <_Z14updatePositionP10particle_td+0x76>
        *(positionDim++) += *(velocityDim++) * timeStep;
   10ab9:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10abd:	48 8d 50 08          	lea    0x8(%rax),%rdx
   10ac1:	48 89 55 f8          	mov    %rdx,-0x8(%rbp)
   10ac5:	f2 0f 10 00          	movsd  (%rax),%xmm0
   10ac9:	66 0f 28 c8          	movapd %xmm0,%xmm1
   10acd:	f2 0f 59 4d d0       	mulsd  -0x30(%rbp),%xmm1
   10ad2:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   10ad6:	48 8d 50 08          	lea    0x8(%rax),%rdx
   10ada:	48 89 55 f0          	mov    %rdx,-0x10(%rbp)
   10ade:	f2 0f 10 00          	movsd  (%rax),%xmm0
   10ae2:	f2 0f 58 c1          	addsd  %xmm1,%xmm0
   10ae6:	f2 0f 11 00          	movsd  %xmm0,(%rax)
    for(int n=0; n<particle->ndim; n++) 
   10aea:	83 45 ec 01          	addl   $0x1,-0x14(%rbp)
   10aee:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10af2:	8b 40 30             	mov    0x30(%rax),%eax
   10af5:	39 45 ec             	cmp    %eax,-0x14(%rbp)
   10af8:	7c bf                	jl     10ab9 <_Z14updatePositionP10particle_td+0x41>
   10afa:	eb 01                	jmp    10afd <_Z14updatePositionP10particle_td+0x85>
        return;
   10afc:	90                   	nop
}
   10afd:	5d                   	pop    %rbp
   10afe:	c3                   	ret    

0000000000010aff <_Z18updateAccelerationP10particle_tS0_>:

void updateAcceleration(particle_t* one, particle_t* another) {
   10aff:	f3 0f 1e fa          	endbr64 
   10b03:	55                   	push   %rbp
   10b04:	48 89 e5             	mov    %rsp,%rbp
   10b07:	48 83 ec 20          	sub    $0x20,%rsp
   10b0b:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   10b0f:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    if (!one->enabled) {
   10b13:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10b17:	0f b6 40 38          	movzbl 0x38(%rax),%eax
   10b1b:	83 f0 01             	xor    $0x1,%eax
   10b1e:	84 c0                	test   %al,%al
   10b20:	74 3a                	je     10b5c <_Z18updateAccelerationP10particle_tS0_+0x5d>
        for(int n=0; n<one->ndim; n++)
   10b22:	c7 45 f8 00 00 00 00 	movl   $0x0,-0x8(%rbp)
   10b29:	eb 20                	jmp    10b4b <_Z18updateAccelerationP10particle_tS0_+0x4c>
            *((one->acceleration)+n) = 0.0;
   10b2b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10b2f:	48 8b 50 18          	mov    0x18(%rax),%rdx
   10b33:	8b 45 f8             	mov    -0x8(%rbp),%eax
   10b36:	48 98                	cltq   
   10b38:	48 c1 e0 03          	shl    $0x3,%rax
   10b3c:	48 01 d0             	add    %rdx,%rax
   10b3f:	66 0f ef c0          	pxor   %xmm0,%xmm0
   10b43:	f2 0f 11 00          	movsd  %xmm0,(%rax)
        for(int n=0; n<one->ndim; n++)
   10b47:	83 45 f8 01          	addl   $0x1,-0x8(%rbp)
   10b4b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10b4f:	8b 40 30             	mov    0x30(%rax),%eax
   10b52:	39 45 f8             	cmp    %eax,-0x8(%rbp)
   10b55:	7c d4                	jl     10b2b <_Z18updateAccelerationP10particle_tS0_+0x2c>
        return;
   10b57:	e9 ed 00 00 00       	jmp    10c49 <_Z18updateAccelerationP10particle_tS0_+0x14a>
    }
    if (!another->enabled) {
   10b5c:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   10b60:	0f b6 40 38          	movzbl 0x38(%rax),%eax
   10b64:	83 f0 01             	xor    $0x1,%eax
   10b67:	84 c0                	test   %al,%al
   10b69:	74 3a                	je     10ba5 <_Z18updateAccelerationP10particle_tS0_+0xa6>
        for(int n=0; n<another->ndim; n++)
   10b6b:	c7 45 fc 00 00 00 00 	movl   $0x0,-0x4(%rbp)
   10b72:	eb 20                	jmp    10b94 <_Z18updateAccelerationP10particle_tS0_+0x95>
            *((another->acceleration)+n) = 0.0;
   10b74:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   10b78:	48 8b 50 18          	mov    0x18(%rax),%rdx
   10b7c:	8b 45 fc             	mov    -0x4(%rbp),%eax
   10b7f:	48 98                	cltq   
   10b81:	48 c1 e0 03          	shl    $0x3,%rax
   10b85:	48 01 d0             	add    %rdx,%rax
   10b88:	66 0f ef c0          	pxor   %xmm0,%xmm0
   10b8c:	f2 0f 11 00          	movsd  %xmm0,(%rax)
        for(int n=0; n<another->ndim; n++)
   10b90:	83 45 fc 01          	addl   $0x1,-0x4(%rbp)
   10b94:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   10b98:	8b 40 30             	mov    0x30(%rax),%eax
   10b9b:	39 45 fc             	cmp    %eax,-0x4(%rbp)
   10b9e:	7c d4                	jl     10b74 <_Z18updateAccelerationP10particle_tS0_+0x75>
        return;
   10ba0:	e9 a4 00 00 00       	jmp    10c49 <_Z18updateAccelerationP10particle_tS0_+0x14a>
    }
    assert(one->updateAcceleration == another->updateAcceleration);
   10ba5:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10ba9:	48 8b 50 28          	mov    0x28(%rax),%rdx
   10bad:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   10bb1:	48 8b 40 28          	mov    0x28(%rax),%rax
   10bb5:	48 39 c2             	cmp    %rax,%rdx
   10bb8:	74 28                	je     10be2 <_Z18updateAccelerationP10particle_tS0_+0xe3>
   10bba:	48 8d 05 a7 5b 00 00 	lea    0x5ba7(%rip),%rax        # 16768 <_ZStL19piecewise_construct+0x8>
   10bc1:	48 89 c1             	mov    %rax,%rcx
   10bc4:	ba 2b 00 00 00       	mov    $0x2b,%edx
   10bc9:	48 8d 05 d0 5b 00 00 	lea    0x5bd0(%rip),%rax        # 167a0 <_ZStL19piecewise_construct+0x40>
   10bd0:	48 89 c6             	mov    %rax,%rsi
   10bd3:	48 8d 05 1e 5c 00 00 	lea    0x5c1e(%rip),%rax        # 167f8 <_ZStL19piecewise_construct+0x98>
   10bda:	48 89 c7             	mov    %rax,%rdi
   10bdd:	e8 6e 95 ff ff       	call   a150 <__assert_fail@plt>
    if (one->updateAcceleration == nullptr) 
   10be2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10be6:	48 8b 40 28          	mov    0x28(%rax),%rax
   10bea:	48 85 c0             	test   %rax,%rax
   10bed:	75 42                	jne    10c31 <_Z18updateAccelerationP10particle_tS0_+0x132>
        std::cerr << "no updateAcceleration registerd on particle: " << one->id << std::endl;
   10bef:	48 8d 05 3a 5c 00 00 	lea    0x5c3a(%rip),%rax        # 16830 <_ZStL19piecewise_construct+0xd0>
   10bf6:	48 89 c6             	mov    %rax,%rsi
   10bf9:	48 8d 05 40 d6 00 00 	lea    0xd640(%rip),%rax        # 1e240 <_ZSt4cerr@GLIBCXX_3.4>
   10c00:	48 89 c7             	mov    %rax,%rdi
   10c03:	e8 c8 97 ff ff       	call   a3d0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
   10c08:	48 89 c2             	mov    %rax,%rdx
   10c0b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10c0f:	48 8b 00             	mov    (%rax),%rax
   10c12:	48 89 c6             	mov    %rax,%rsi
   10c15:	48 89 d7             	mov    %rdx,%rdi
   10c18:	e8 d3 92 ff ff       	call   9ef0 <_ZNSolsEm@plt>
   10c1d:	48 8b 15 ac c3 00 00 	mov    0xc3ac(%rip),%rdx        # 1cfd0 <_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_@GLIBCXX_3.4>
   10c24:	48 89 d6             	mov    %rdx,%rsi
   10c27:	48 89 c7             	mov    %rax,%rdi
   10c2a:	e8 e1 97 ff ff       	call   a410 <_ZNSolsEPFRSoS_E@plt>
   10c2f:	eb 18                	jmp    10c49 <_Z18updateAccelerationP10particle_tS0_+0x14a>
    else
        (*(one->updateAcceleration))(one, another);
   10c31:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10c35:	48 8b 48 28          	mov    0x28(%rax),%rcx
   10c39:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
   10c3d:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10c41:	48 89 d6             	mov    %rdx,%rsi
   10c44:	48 89 c7             	mov    %rax,%rdi
   10c47:	ff d1                	call   *%rcx
}
   10c49:	c9                   	leave  
   10c4a:	c3                   	ret    

0000000000010c4b <_Z28registerUpdateAccelerationFnP10particle_tPFvS0_S0_E>:

void registerUpdateAccelerationFn(particle_t* one, void (*fn)(particle_t*, particle_t*)) {
   10c4b:	f3 0f 1e fa          	endbr64 
   10c4f:	55                   	push   %rbp
   10c50:	48 89 e5             	mov    %rsp,%rbp
   10c53:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10c57:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    one->updateAcceleration = fn;
   10c5b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10c5f:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   10c63:	48 89 50 28          	mov    %rdx,0x28(%rax)
}
   10c67:	90                   	nop
   10c68:	5d                   	pop    %rbp
   10c69:	c3                   	ret    

0000000000010c6a <_Z29registerUpdateAccelerationFnsP10particle_tiPFvS0_S0_E>:

void registerUpdateAccelerationFns(particle_t* particles, int nParticles, void (*fn)(particle_t*, particle_t*)) {
   10c6a:	f3 0f 1e fa          	endbr64 
   10c6e:	55                   	push   %rbp
   10c6f:	48 89 e5             	mov    %rsp,%rbp
   10c72:	48 83 ec 28          	sub    $0x28,%rsp
   10c76:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   10c7a:	89 75 e4             	mov    %esi,-0x1c(%rbp)
   10c7d:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    particle_t* pt = particles;
   10c81:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10c85:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    for (int n=0; n<nParticles; n++)
   10c89:	c7 45 f4 00 00 00 00 	movl   $0x0,-0xc(%rbp)
   10c90:	eb 1f                	jmp    10cb1 <_Z29registerUpdateAccelerationFnsP10particle_tiPFvS0_S0_E+0x47>
        registerUpdateAccelerationFn(pt++, fn);
   10c92:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10c96:	48 8d 50 40          	lea    0x40(%rax),%rdx
   10c9a:	48 89 55 f8          	mov    %rdx,-0x8(%rbp)
   10c9e:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
   10ca2:	48 89 d6             	mov    %rdx,%rsi
   10ca5:	48 89 c7             	mov    %rax,%rdi
   10ca8:	e8 9e ff ff ff       	call   10c4b <_Z28registerUpdateAccelerationFnP10particle_tPFvS0_S0_E>
    for (int n=0; n<nParticles; n++)
   10cad:	83 45 f4 01          	addl   $0x1,-0xc(%rbp)
   10cb1:	8b 45 f4             	mov    -0xc(%rbp),%eax
   10cb4:	3b 45 e4             	cmp    -0x1c(%rbp),%eax
   10cb7:	7c d9                	jl     10c92 <_Z29registerUpdateAccelerationFnsP10particle_tiPFvS0_S0_E+0x28>
}
   10cb9:	90                   	nop
   10cba:	90                   	nop
   10cbb:	c9                   	leave  
   10cbc:	c3                   	ret    

0000000000010cbd <_Z29registerUpdateAccelerationFnsP17chunk_particles_tPFvP10particle_tS2_E>:

void registerUpdateAccelerationFns(chunk_particles_t* particleChunk, void (*fn)(particle_t*, particle_t*)) {
   10cbd:	f3 0f 1e fa          	endbr64 
   10cc1:	55                   	push   %rbp
   10cc2:	48 89 e5             	mov    %rsp,%rbp
   10cc5:	48 83 ec 10          	sub    $0x10,%rsp
   10cc9:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10ccd:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    registerUpdateAccelerationFns(particleChunk->particles, particleChunk->nParticle, fn);
   10cd1:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10cd5:	8b 48 08             	mov    0x8(%rax),%ecx
   10cd8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10cdc:	48 8b 00             	mov    (%rax),%rax
   10cdf:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   10ce3:	89 ce                	mov    %ecx,%esi
   10ce5:	48 89 c7             	mov    %rax,%rdi
   10ce8:	e8 7d ff ff ff       	call   10c6a <_Z29registerUpdateAccelerationFnsP10particle_tiPFvS0_S0_E>
}
   10ced:	90                   	nop
   10cee:	c9                   	leave  
   10cef:	c3                   	ret    

0000000000010cf0 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi>:


void load_particles(const std::string filename, particle_t** particles, int* nParticles) {
   10cf0:	f3 0f 1e fa          	endbr64 
   10cf4:	55                   	push   %rbp
   10cf5:	48 89 e5             	mov    %rsp,%rbp
   10cf8:	53                   	push   %rbx
   10cf9:	48 81 ec 48 04 00 00 	sub    $0x448,%rsp
   10d00:	48 89 bd c8 fb ff ff 	mov    %rdi,-0x438(%rbp)
   10d07:	48 89 b5 c0 fb ff ff 	mov    %rsi,-0x440(%rbp)
   10d0e:	48 89 95 b8 fb ff ff 	mov    %rdx,-0x448(%rbp)
   10d15:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   10d1c:	00 00 
   10d1e:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   10d22:	31 c0                	xor    %eax,%eax
    std::ifstream file(filename, std::ios::in);
   10d24:	48 8b 8d c8 fb ff ff 	mov    -0x438(%rbp),%rcx
   10d2b:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   10d32:	ba 08 00 00 00       	mov    $0x8,%edx
   10d37:	48 89 ce             	mov    %rcx,%rsi
   10d3a:	48 89 c7             	mov    %rax,%rdi
   10d3d:	e8 5e 9b ff ff       	call   a8a0 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEEC1ERKNSt7__cxx1112basic_stringIcS1_SaIcEEESt13_Ios_Openmode@plt>
    if(!file.is_open()) {
   10d42:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   10d49:	48 89 c7             	mov    %rax,%rdi
   10d4c:	e8 1f 9b ff ff       	call   a870 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEE7is_openEv@plt>
   10d51:	83 f0 01             	xor    $0x1,%eax
   10d54:	84 c0                	test   %al,%al
   10d56:	74 63                	je     10dbb <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0xcb>
        std::fprintf(stderr, "[%s:%d]: can not open the file %s\n", __FILE__, __LINE__, filename);
   10d58:	48 8b 95 c8 fb ff ff 	mov    -0x438(%rbp),%rdx
   10d5f:	48 8d 85 50 fc ff ff 	lea    -0x3b0(%rbp),%rax
   10d66:	48 89 d6             	mov    %rdx,%rsi
   10d69:	48 89 c7             	mov    %rax,%rdi
   10d6c:	e8 6f 92 ff ff       	call   9fe0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
   10d71:	48 8b 05 a8 d4 00 00 	mov    0xd4a8(%rip),%rax        # 1e220 <stderr@GLIBC_2.2.5>
   10d78:	48 8d 95 50 fc ff ff 	lea    -0x3b0(%rbp),%rdx
   10d7f:	49 89 d0             	mov    %rdx,%r8
   10d82:	b9 44 00 00 00       	mov    $0x44,%ecx
   10d87:	48 8d 15 12 5a 00 00 	lea    0x5a12(%rip),%rdx        # 167a0 <_ZStL19piecewise_construct+0x40>
   10d8e:	48 8d 35 cb 5a 00 00 	lea    0x5acb(%rip),%rsi        # 16860 <_ZStL19piecewise_construct+0x100>
   10d95:	48 89 c7             	mov    %rax,%rdi
   10d98:	b8 00 00 00 00       	mov    $0x0,%eax
   10d9d:	e8 be 99 ff ff       	call   a760 <fprintf@plt>
   10da2:	48 8d 85 50 fc ff ff 	lea    -0x3b0(%rbp),%rax
   10da9:	48 89 c7             	mov    %rax,%rdi
   10dac:	e8 6f 93 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
        exit(EXIT_FAILURE);
   10db1:	bf 01 00 00 00       	mov    $0x1,%edi
   10db6:	e8 85 97 ff ff       	call   a540 <exit@plt>
    }
    int space_dimension = -1;
   10dbb:	c7 85 d0 fb ff ff ff 	movl   $0xffffffff,-0x430(%rbp)
   10dc2:	ff ff ff 
    std::string line;
   10dc5:	48 8d 85 30 fc ff ff 	lea    -0x3d0(%rbp),%rax
   10dcc:	48 89 c7             	mov    %rax,%rdi
   10dcf:	e8 cc 98 ff ff       	call   a6a0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1Ev@plt>
    std::getline(file, line);
   10dd4:	48 8d 95 30 fc ff ff 	lea    -0x3d0(%rbp),%rdx
   10ddb:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   10de2:	48 89 d6             	mov    %rdx,%rsi
   10de5:	48 89 c7             	mov    %rax,%rdi
   10de8:	e8 93 90 ff ff       	call   9e80 <_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE@plt>
    std::sscanf(line.c_str(), "%d %d", nParticles, &space_dimension);
   10ded:	48 8d 85 30 fc ff ff 	lea    -0x3d0(%rbp),%rax
   10df4:	48 89 c7             	mov    %rax,%rdi
   10df7:	e8 c4 91 ff ff       	call   9fc0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE5c_strEv@plt>
   10dfc:	48 89 c7             	mov    %rax,%rdi
   10dff:	48 8d 95 d0 fb ff ff 	lea    -0x430(%rbp),%rdx
   10e06:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   10e0d:	48 89 d1             	mov    %rdx,%rcx
   10e10:	48 89 c2             	mov    %rax,%rdx
   10e13:	48 8d 05 69 5a 00 00 	lea    0x5a69(%rip),%rax        # 16883 <_ZStL19piecewise_construct+0x123>
   10e1a:	48 89 c6             	mov    %rax,%rsi
   10e1d:	b8 00 00 00 00       	mov    $0x0,%eax
   10e22:	e8 99 95 ff ff       	call   a3c0 <__isoc99_sscanf@plt>
    particle_t* particles_ = static_cast<particle_t*>(malloc(*nParticles * sizeof(particle_t)));
   10e27:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   10e2e:	8b 00                	mov    (%rax),%eax
   10e30:	48 98                	cltq   
   10e32:	48 c1 e0 06          	shl    $0x6,%rax
   10e36:	48 89 c7             	mov    %rax,%rdi
   10e39:	e8 f2 97 ff ff       	call   a630 <malloc@plt>
   10e3e:	48 89 85 00 fc ff ff 	mov    %rax,-0x400(%rbp)
    double* pos = static_cast<double*>(malloc((*nParticles) * space_dimension * sizeof(double)));
   10e45:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   10e4c:	8b 10                	mov    (%rax),%edx
   10e4e:	8b 85 d0 fb ff ff    	mov    -0x430(%rbp),%eax
   10e54:	0f af c2             	imul   %edx,%eax
   10e57:	48 98                	cltq   
   10e59:	48 c1 e0 03          	shl    $0x3,%rax
   10e5d:	48 89 c7             	mov    %rax,%rdi
   10e60:	e8 cb 97 ff ff       	call   a630 <malloc@plt>
   10e65:	48 89 85 08 fc ff ff 	mov    %rax,-0x3f8(%rbp)
    assert(pos != nullptr);
   10e6c:	48 83 bd 08 fc ff ff 	cmpq   $0x0,-0x3f8(%rbp)
   10e73:	00 
   10e74:	75 28                	jne    10e9e <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x1ae>
   10e76:	48 8d 05 13 5a 00 00 	lea    0x5a13(%rip),%rax        # 16890 <_ZStL19piecewise_construct+0x130>
   10e7d:	48 89 c1             	mov    %rax,%rcx
   10e80:	ba 4d 00 00 00       	mov    $0x4d,%edx
   10e85:	48 8d 05 14 59 00 00 	lea    0x5914(%rip),%rax        # 167a0 <_ZStL19piecewise_construct+0x40>
   10e8c:	48 89 c6             	mov    %rax,%rsi
   10e8f:	48 8d 05 2f 5a 00 00 	lea    0x5a2f(%rip),%rax        # 168c5 <_ZStL19piecewise_construct+0x165>
   10e96:	48 89 c7             	mov    %rax,%rdi
   10e99:	e8 b2 92 ff ff       	call   a150 <__assert_fail@plt>
    double* vel = static_cast<double*>(malloc((*nParticles) * space_dimension * sizeof(double)));
   10e9e:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   10ea5:	8b 10                	mov    (%rax),%edx
   10ea7:	8b 85 d0 fb ff ff    	mov    -0x430(%rbp),%eax
   10ead:	0f af c2             	imul   %edx,%eax
   10eb0:	48 98                	cltq   
   10eb2:	48 c1 e0 03          	shl    $0x3,%rax
   10eb6:	48 89 c7             	mov    %rax,%rdi
   10eb9:	e8 72 97 ff ff       	call   a630 <malloc@plt>
   10ebe:	48 89 85 10 fc ff ff 	mov    %rax,-0x3f0(%rbp)
    assert(vel != nullptr);
   10ec5:	48 83 bd 10 fc ff ff 	cmpq   $0x0,-0x3f0(%rbp)
   10ecc:	00 
   10ecd:	75 28                	jne    10ef7 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x207>
   10ecf:	48 8d 05 ba 59 00 00 	lea    0x59ba(%rip),%rax        # 16890 <_ZStL19piecewise_construct+0x130>
   10ed6:	48 89 c1             	mov    %rax,%rcx
   10ed9:	ba 4f 00 00 00       	mov    $0x4f,%edx
   10ede:	48 8d 05 bb 58 00 00 	lea    0x58bb(%rip),%rax        # 167a0 <_ZStL19piecewise_construct+0x40>
   10ee5:	48 89 c6             	mov    %rax,%rsi
   10ee8:	48 8d 05 e5 59 00 00 	lea    0x59e5(%rip),%rax        # 168d4 <_ZStL19piecewise_construct+0x174>
   10eef:	48 89 c7             	mov    %rax,%rdi
   10ef2:	e8 59 92 ff ff       	call   a150 <__assert_fail@plt>
    double* acc = static_cast<double*>(malloc((*nParticles) * space_dimension * sizeof(double)));
   10ef7:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   10efe:	8b 10                	mov    (%rax),%edx
   10f00:	8b 85 d0 fb ff ff    	mov    -0x430(%rbp),%eax
   10f06:	0f af c2             	imul   %edx,%eax
   10f09:	48 98                	cltq   
   10f0b:	48 c1 e0 03          	shl    $0x3,%rax
   10f0f:	48 89 c7             	mov    %rax,%rdi
   10f12:	e8 19 97 ff ff       	call   a630 <malloc@plt>
   10f17:	48 89 85 18 fc ff ff 	mov    %rax,-0x3e8(%rbp)
    assert(acc != nullptr);
   10f1e:	48 83 bd 18 fc ff ff 	cmpq   $0x0,-0x3e8(%rbp)
   10f25:	00 
   10f26:	75 28                	jne    10f50 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x260>
   10f28:	48 8d 05 61 59 00 00 	lea    0x5961(%rip),%rax        # 16890 <_ZStL19piecewise_construct+0x130>
   10f2f:	48 89 c1             	mov    %rax,%rcx
   10f32:	ba 51 00 00 00       	mov    $0x51,%edx
   10f37:	48 8d 05 62 58 00 00 	lea    0x5862(%rip),%rax        # 167a0 <_ZStL19piecewise_construct+0x40>
   10f3e:	48 89 c6             	mov    %rax,%rsi
   10f41:	48 8d 05 9b 59 00 00 	lea    0x599b(%rip),%rax        # 168e3 <_ZStL19piecewise_construct+0x183>
   10f48:	48 89 c7             	mov    %rax,%rdi
   10f4b:	e8 00 92 ff ff       	call   a150 <__assert_fail@plt>
    double* features = nullptr;
   10f50:	48 c7 85 f0 fb ff ff 	movq   $0x0,-0x410(%rbp)
   10f57:	00 00 00 00 
    particle_t* pp = particles_;
   10f5b:	48 8b 85 00 fc ff ff 	mov    -0x400(%rbp),%rax
   10f62:	48 89 85 f8 fb ff ff 	mov    %rax,-0x408(%rbp)
    for(int n=0; n<(*nParticles); n++) {
   10f69:	c7 85 d4 fb ff ff 00 	movl   $0x0,-0x42c(%rbp)
   10f70:	00 00 00 
   10f73:	e9 da 00 00 00       	jmp    11052 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x362>
        pp->position = pos+n*space_dimension;
   10f78:	8b 85 d0 fb ff ff    	mov    -0x430(%rbp),%eax
   10f7e:	0f af 85 d4 fb ff ff 	imul   -0x42c(%rbp),%eax
   10f85:	48 98                	cltq   
   10f87:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   10f8e:	00 
   10f8f:	48 8b 85 08 fc ff ff 	mov    -0x3f8(%rbp),%rax
   10f96:	48 01 c2             	add    %rax,%rdx
   10f99:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   10fa0:	48 89 50 08          	mov    %rdx,0x8(%rax)
        pp->velocity = vel+n*space_dimension;
   10fa4:	8b 85 d0 fb ff ff    	mov    -0x430(%rbp),%eax
   10faa:	0f af 85 d4 fb ff ff 	imul   -0x42c(%rbp),%eax
   10fb1:	48 98                	cltq   
   10fb3:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   10fba:	00 
   10fbb:	48 8b 85 10 fc ff ff 	mov    -0x3f0(%rbp),%rax
   10fc2:	48 01 c2             	add    %rax,%rdx
   10fc5:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   10fcc:	48 89 50 10          	mov    %rdx,0x10(%rax)
        pp->acceleration = acc+n*space_dimension;
   10fd0:	8b 85 d0 fb ff ff    	mov    -0x430(%rbp),%eax
   10fd6:	0f af 85 d4 fb ff ff 	imul   -0x42c(%rbp),%eax
   10fdd:	48 98                	cltq   
   10fdf:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   10fe6:	00 
   10fe7:	48 8b 85 18 fc ff ff 	mov    -0x3e8(%rbp),%rax
   10fee:	48 01 c2             	add    %rax,%rdx
   10ff1:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   10ff8:	48 89 50 18          	mov    %rdx,0x18(%rax)
        pp->ndim = space_dimension;
   10ffc:	8b 95 d0 fb ff ff    	mov    -0x430(%rbp),%edx
   11002:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   11009:	89 50 30             	mov    %edx,0x30(%rax)
        pp->updateAcceleration = nullptr;
   1100c:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   11013:	48 c7 40 28 00 00 00 	movq   $0x0,0x28(%rax)
   1101a:	00 
        pp->features = nullptr;
   1101b:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   11022:	48 c7 40 20 00 00 00 	movq   $0x0,0x20(%rax)
   11029:	00 
        pp->nfeat = -1;
   1102a:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   11031:	c7 40 34 ff ff ff ff 	movl   $0xffffffff,0x34(%rax)
        pp->enabled=true;
   11038:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   1103f:	c6 40 38 01          	movb   $0x1,0x38(%rax)
        pp++;
   11043:	48 83 85 f8 fb ff ff 	addq   $0x40,-0x408(%rbp)
   1104a:	40 
    for(int n=0; n<(*nParticles); n++) {
   1104b:	83 85 d4 fb ff ff 01 	addl   $0x1,-0x42c(%rbp)
   11052:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   11059:	8b 00                	mov    (%rax),%eax
   1105b:	39 85 d4 fb ff ff    	cmp    %eax,-0x42c(%rbp)
   11061:	0f 8c 11 ff ff ff    	jl     10f78 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x288>
    }

    pp = particles_;
   11067:	48 8b 85 00 fc ff ff 	mov    -0x400(%rbp),%rax
   1106e:	48 89 85 f8 fb ff ff 	mov    %rax,-0x408(%rbp)
    std::stringstream lineStream;
   11075:	48 8d 85 50 fc ff ff 	lea    -0x3b0(%rbp),%rax
   1107c:	48 89 c7             	mov    %rax,%rdi
   1107f:	e8 4c 91 ff ff       	call   a1d0 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEC1Ev@plt>
    double foo;
    int n = 0;
   11084:	c7 85 d8 fb ff ff 00 	movl   $0x0,-0x428(%rbp)
   1108b:	00 00 00 
    int nfeat = -1;
   1108e:	c7 85 dc fb ff ff ff 	movl   $0xffffffff,-0x424(%rbp)
   11095:	ff ff ff 
    while(std::getline(file, line) && n<(*nParticles)) {
   11098:	e9 60 03 00 00       	jmp    113fd <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x70d>
        char* bar;
        bar = strtok(const_cast<char*>(line.c_str()), " ");
   1109d:	48 8d 85 30 fc ff ff 	lea    -0x3d0(%rbp),%rax
   110a4:	48 89 c7             	mov    %rax,%rdi
   110a7:	e8 14 8f ff ff       	call   9fc0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE5c_strEv@plt>
   110ac:	48 89 c2             	mov    %rax,%rdx
   110af:	48 8d 05 3c 58 00 00 	lea    0x583c(%rip),%rax        # 168f2 <_ZStL19piecewise_construct+0x192>
   110b6:	48 89 c6             	mov    %rax,%rsi
   110b9:	48 89 d7             	mov    %rdx,%rdi
   110bc:	e8 6f 97 ff ff       	call   a830 <strtok@plt>
   110c1:	48 89 85 20 fc ff ff 	mov    %rax,-0x3e0(%rbp)
        foo = atof(bar);
   110c8:	48 8b 85 20 fc ff ff 	mov    -0x3e0(%rbp),%rax
   110cf:	48 89 c7             	mov    %rax,%rdi
   110d2:	e8 a9 94 ff ff       	call   a580 <atof@plt>
   110d7:	66 48 0f 7e c0       	movq   %xmm0,%rax
   110dc:	48 89 85 28 fc ff ff 	mov    %rax,-0x3d8(%rbp)
        pp->id = static_cast<uint64_t>(foo);
   110e3:	f2 0f 10 85 28 fc ff 	movsd  -0x3d8(%rbp),%xmm0
   110ea:	ff 
   110eb:	66 0f 2f 05 65 59 00 	comisd 0x5965(%rip),%xmm0        # 16a58 <_ZStL19piecewise_construct+0x2f8>
   110f2:	00 
   110f3:	73 0f                	jae    11104 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x414>
   110f5:	f2 0f 10 85 28 fc ff 	movsd  -0x3d8(%rbp),%xmm0
   110fc:	ff 
   110fd:	f2 48 0f 2c c0       	cvttsd2si %xmm0,%rax
   11102:	eb 26                	jmp    1112a <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x43a>
   11104:	f2 0f 10 85 28 fc ff 	movsd  -0x3d8(%rbp),%xmm0
   1110b:	ff 
   1110c:	f2 0f 10 0d 44 59 00 	movsd  0x5944(%rip),%xmm1        # 16a58 <_ZStL19piecewise_construct+0x2f8>
   11113:	00 
   11114:	f2 0f 5c c1          	subsd  %xmm1,%xmm0
   11118:	f2 48 0f 2c c0       	cvttsd2si %xmm0,%rax
   1111d:	48 ba 00 00 00 00 00 	movabs $0x8000000000000000,%rdx
   11124:	00 00 80 
   11127:	48 31 d0             	xor    %rdx,%rax
   1112a:	48 8b 95 f8 fb ff ff 	mov    -0x408(%rbp),%rdx
   11131:	48 89 02             	mov    %rax,(%rdx)
        for(int nn=0; nn<(pp->ndim); nn++) {
   11134:	c7 85 e0 fb ff ff 00 	movl   $0x0,-0x420(%rbp)
   1113b:	00 00 00 
   1113e:	eb 54                	jmp    11194 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x4a4>
            bar = strtok(NULL, " ");
   11140:	48 8d 05 ab 57 00 00 	lea    0x57ab(%rip),%rax        # 168f2 <_ZStL19piecewise_construct+0x192>
   11147:	48 89 c6             	mov    %rax,%rsi
   1114a:	bf 00 00 00 00       	mov    $0x0,%edi
   1114f:	e8 dc 96 ff ff       	call   a830 <strtok@plt>
   11154:	48 89 85 20 fc ff ff 	mov    %rax,-0x3e0(%rbp)
            *(pp->position + nn) = atof(bar);
   1115b:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   11162:	48 8b 50 08          	mov    0x8(%rax),%rdx
   11166:	8b 85 e0 fb ff ff    	mov    -0x420(%rbp),%eax
   1116c:	48 98                	cltq   
   1116e:	48 c1 e0 03          	shl    $0x3,%rax
   11172:	48 8d 1c 02          	lea    (%rdx,%rax,1),%rbx
   11176:	48 8b 85 20 fc ff ff 	mov    -0x3e0(%rbp),%rax
   1117d:	48 89 c7             	mov    %rax,%rdi
   11180:	e8 fb 93 ff ff       	call   a580 <atof@plt>
   11185:	66 48 0f 7e c0       	movq   %xmm0,%rax
   1118a:	48 89 03             	mov    %rax,(%rbx)
        for(int nn=0; nn<(pp->ndim); nn++) {
   1118d:	83 85 e0 fb ff ff 01 	addl   $0x1,-0x420(%rbp)
   11194:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   1119b:	8b 40 30             	mov    0x30(%rax),%eax
   1119e:	39 85 e0 fb ff ff    	cmp    %eax,-0x420(%rbp)
   111a4:	7c 9a                	jl     11140 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x450>
        }
        for(int nn=0; nn<(pp->ndim); nn++) {
   111a6:	c7 85 e4 fb ff ff 00 	movl   $0x0,-0x41c(%rbp)
   111ad:	00 00 00 
   111b0:	eb 54                	jmp    11206 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x516>
            bar = strtok(NULL, " ");
   111b2:	48 8d 05 39 57 00 00 	lea    0x5739(%rip),%rax        # 168f2 <_ZStL19piecewise_construct+0x192>
   111b9:	48 89 c6             	mov    %rax,%rsi
   111bc:	bf 00 00 00 00       	mov    $0x0,%edi
   111c1:	e8 6a 96 ff ff       	call   a830 <strtok@plt>
   111c6:	48 89 85 20 fc ff ff 	mov    %rax,-0x3e0(%rbp)
            *(pp->velocity + nn) = atof(bar);
   111cd:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   111d4:	48 8b 50 10          	mov    0x10(%rax),%rdx
   111d8:	8b 85 e4 fb ff ff    	mov    -0x41c(%rbp),%eax
   111de:	48 98                	cltq   
   111e0:	48 c1 e0 03          	shl    $0x3,%rax
   111e4:	48 8d 1c 02          	lea    (%rdx,%rax,1),%rbx
   111e8:	48 8b 85 20 fc ff ff 	mov    -0x3e0(%rbp),%rax
   111ef:	48 89 c7             	mov    %rax,%rdi
   111f2:	e8 89 93 ff ff       	call   a580 <atof@plt>
   111f7:	66 48 0f 7e c0       	movq   %xmm0,%rax
   111fc:	48 89 03             	mov    %rax,(%rbx)
        for(int nn=0; nn<(pp->ndim); nn++) {
   111ff:	83 85 e4 fb ff ff 01 	addl   $0x1,-0x41c(%rbp)
   11206:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   1120d:	8b 40 30             	mov    0x30(%rax),%eax
   11210:	39 85 e4 fb ff ff    	cmp    %eax,-0x41c(%rbp)
   11216:	7c 9a                	jl     111b2 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x4c2>
        }
        for(int nn=0; nn<(pp->ndim); nn++) {
   11218:	c7 85 e8 fb ff ff 00 	movl   $0x0,-0x418(%rbp)
   1121f:	00 00 00 
   11222:	eb 54                	jmp    11278 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x588>
            bar = strtok(NULL, " ");
   11224:	48 8d 05 c7 56 00 00 	lea    0x56c7(%rip),%rax        # 168f2 <_ZStL19piecewise_construct+0x192>
   1122b:	48 89 c6             	mov    %rax,%rsi
   1122e:	bf 00 00 00 00       	mov    $0x0,%edi
   11233:	e8 f8 95 ff ff       	call   a830 <strtok@plt>
   11238:	48 89 85 20 fc ff ff 	mov    %rax,-0x3e0(%rbp)
            *(pp->acceleration + nn) = atof(bar);
   1123f:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   11246:	48 8b 50 18          	mov    0x18(%rax),%rdx
   1124a:	8b 85 e8 fb ff ff    	mov    -0x418(%rbp),%eax
   11250:	48 98                	cltq   
   11252:	48 c1 e0 03          	shl    $0x3,%rax
   11256:	48 8d 1c 02          	lea    (%rdx,%rax,1),%rbx
   1125a:	48 8b 85 20 fc ff ff 	mov    -0x3e0(%rbp),%rax
   11261:	48 89 c7             	mov    %rax,%rdi
   11264:	e8 17 93 ff ff       	call   a580 <atof@plt>
   11269:	66 48 0f 7e c0       	movq   %xmm0,%rax
   1126e:	48 89 03             	mov    %rax,(%rbx)
        for(int nn=0; nn<(pp->ndim); nn++) {
   11271:	83 85 e8 fb ff ff 01 	addl   $0x1,-0x418(%rbp)
   11278:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   1127f:	8b 40 30             	mov    0x30(%rax),%eax
   11282:	39 85 e8 fb ff ff    	cmp    %eax,-0x418(%rbp)
   11288:	7c 9a                	jl     11224 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x534>
        }
        bar = strtok(NULL, " ");
   1128a:	48 8d 05 61 56 00 00 	lea    0x5661(%rip),%rax        # 168f2 <_ZStL19piecewise_construct+0x192>
   11291:	48 89 c6             	mov    %rax,%rsi
   11294:	bf 00 00 00 00       	mov    $0x0,%edi
   11299:	e8 92 95 ff ff       	call   a830 <strtok@plt>
   1129e:	48 89 85 20 fc ff ff 	mov    %rax,-0x3e0(%rbp)
        foo = atof(bar);
   112a5:	48 8b 85 20 fc ff ff 	mov    -0x3e0(%rbp),%rax
   112ac:	48 89 c7             	mov    %rax,%rdi
   112af:	e8 cc 92 ff ff       	call   a580 <atof@plt>
   112b4:	66 48 0f 7e c0       	movq   %xmm0,%rax
   112b9:	48 89 85 28 fc ff ff 	mov    %rax,-0x3d8(%rbp)
        pp->nfeat = static_cast<int>(foo);
   112c0:	f2 0f 10 85 28 fc ff 	movsd  -0x3d8(%rbp),%xmm0
   112c7:	ff 
   112c8:	f2 0f 2c d0          	cvttsd2si %xmm0,%edx
   112cc:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   112d3:	89 50 34             	mov    %edx,0x34(%rax)
        if (nfeat == -1) {
   112d6:	83 bd dc fb ff ff ff 	cmpl   $0xffffffff,-0x424(%rbp)
   112dd:	75 37                	jne    11316 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x626>
            nfeat = pp->nfeat;
   112df:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   112e6:	8b 40 34             	mov    0x34(%rax),%eax
   112e9:	89 85 dc fb ff ff    	mov    %eax,-0x424(%rbp)
            features = static_cast<double*>(malloc(nfeat*(*nParticles)*sizeof(double)));
   112ef:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   112f6:	8b 00                	mov    (%rax),%eax
   112f8:	0f af 85 dc fb ff ff 	imul   -0x424(%rbp),%eax
   112ff:	48 98                	cltq   
   11301:	48 c1 e0 03          	shl    $0x3,%rax
   11305:	48 89 c7             	mov    %rax,%rdi
   11308:	e8 23 93 ff ff       	call   a630 <malloc@plt>
   1130d:	48 89 85 f0 fb ff ff 	mov    %rax,-0x410(%rbp)
   11314:	eb 3a                	jmp    11350 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x660>
        }
        else {
            assert(nfeat==(pp->nfeat));
   11316:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   1131d:	8b 40 34             	mov    0x34(%rax),%eax
   11320:	39 85 dc fb ff ff    	cmp    %eax,-0x424(%rbp)
   11326:	74 28                	je     11350 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x660>
   11328:	48 8d 05 61 55 00 00 	lea    0x5561(%rip),%rax        # 16890 <_ZStL19piecewise_construct+0x130>
   1132f:	48 89 c1             	mov    %rax,%rcx
   11332:	ba 7e 00 00 00       	mov    $0x7e,%edx
   11337:	48 8d 05 62 54 00 00 	lea    0x5462(%rip),%rax        # 167a0 <_ZStL19piecewise_construct+0x40>
   1133e:	48 89 c6             	mov    %rax,%rsi
   11341:	48 8d 05 ac 55 00 00 	lea    0x55ac(%rip),%rax        # 168f4 <_ZStL19piecewise_construct+0x194>
   11348:	48 89 c7             	mov    %rax,%rdi
   1134b:	e8 00 8e ff ff       	call   a150 <__assert_fail@plt>
        }
        pp->features = features + n*nfeat;
   11350:	8b 85 d8 fb ff ff    	mov    -0x428(%rbp),%eax
   11356:	0f af 85 dc fb ff ff 	imul   -0x424(%rbp),%eax
   1135d:	48 98                	cltq   
   1135f:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   11366:	00 
   11367:	48 8b 85 f0 fb ff ff 	mov    -0x410(%rbp),%rax
   1136e:	48 01 c2             	add    %rax,%rdx
   11371:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   11378:	48 89 50 20          	mov    %rdx,0x20(%rax)
        for(int nn=0; nn<(pp->nfeat); nn++) {
   1137c:	c7 85 ec fb ff ff 00 	movl   $0x0,-0x414(%rbp)
   11383:	00 00 00 
   11386:	eb 54                	jmp    113dc <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x6ec>
            bar = strtok(NULL, " ");
   11388:	48 8d 05 63 55 00 00 	lea    0x5563(%rip),%rax        # 168f2 <_ZStL19piecewise_construct+0x192>
   1138f:	48 89 c6             	mov    %rax,%rsi
   11392:	bf 00 00 00 00       	mov    $0x0,%edi
   11397:	e8 94 94 ff ff       	call   a830 <strtok@plt>
   1139c:	48 89 85 20 fc ff ff 	mov    %rax,-0x3e0(%rbp)
            *(pp->features + nn) = atof(bar);
   113a3:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   113aa:	48 8b 50 20          	mov    0x20(%rax),%rdx
   113ae:	8b 85 ec fb ff ff    	mov    -0x414(%rbp),%eax
   113b4:	48 98                	cltq   
   113b6:	48 c1 e0 03          	shl    $0x3,%rax
   113ba:	48 8d 1c 02          	lea    (%rdx,%rax,1),%rbx
   113be:	48 8b 85 20 fc ff ff 	mov    -0x3e0(%rbp),%rax
   113c5:	48 89 c7             	mov    %rax,%rdi
   113c8:	e8 b3 91 ff ff       	call   a580 <atof@plt>
   113cd:	66 48 0f 7e c0       	movq   %xmm0,%rax
   113d2:	48 89 03             	mov    %rax,(%rbx)
        for(int nn=0; nn<(pp->nfeat); nn++) {
   113d5:	83 85 ec fb ff ff 01 	addl   $0x1,-0x414(%rbp)
   113dc:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   113e3:	8b 40 34             	mov    0x34(%rax),%eax
   113e6:	39 85 ec fb ff ff    	cmp    %eax,-0x414(%rbp)
   113ec:	7c 9a                	jl     11388 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x698>
        }
        n++;
   113ee:	83 85 d8 fb ff ff 01 	addl   $0x1,-0x428(%rbp)
        pp++;
   113f5:	48 83 85 f8 fb ff ff 	addq   $0x40,-0x408(%rbp)
   113fc:	40 
    while(std::getline(file, line) && n<(*nParticles)) {
   113fd:	48 8d 95 30 fc ff ff 	lea    -0x3d0(%rbp),%rdx
   11404:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   1140b:	48 89 d6             	mov    %rdx,%rsi
   1140e:	48 89 c7             	mov    %rax,%rdi
   11411:	e8 6a 8a ff ff       	call   9e80 <_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE@plt>
   11416:	48 8b 10             	mov    (%rax),%rdx
   11419:	48 83 ea 18          	sub    $0x18,%rdx
   1141d:	48 8b 12             	mov    (%rdx),%rdx
   11420:	48 01 d0             	add    %rdx,%rax
   11423:	48 89 c7             	mov    %rax,%rdi
   11426:	e8 75 91 ff ff       	call   a5a0 <_ZNKSt9basic_iosIcSt11char_traitsIcEEcvbEv@plt>
   1142b:	84 c0                	test   %al,%al
   1142d:	74 18                	je     11447 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x757>
   1142f:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   11436:	8b 00                	mov    (%rax),%eax
   11438:	39 85 d8 fb ff ff    	cmp    %eax,-0x428(%rbp)
   1143e:	7d 07                	jge    11447 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x757>
   11440:	b8 01 00 00 00       	mov    $0x1,%eax
   11445:	eb 05                	jmp    1144c <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x75c>
   11447:	b8 00 00 00 00       	mov    $0x0,%eax
   1144c:	84 c0                	test   %al,%al
   1144e:	0f 85 49 fc ff ff    	jne    1109d <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x3ad>
    }
    file.close();
   11454:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   1145b:	48 89 c7             	mov    %rax,%rdi
   1145e:	e8 4d 8a ff ff       	call   9eb0 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEE5closeEv@plt>
    *particles = particles_;
   11463:	48 8b 85 c0 fb ff ff 	mov    -0x440(%rbp),%rax
   1146a:	48 8b 95 00 fc ff ff 	mov    -0x400(%rbp),%rdx
   11471:	48 89 10             	mov    %rdx,(%rax)
    if(n != *nParticles) {
   11474:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   1147b:	8b 00                	mov    (%rax),%eax
   1147d:	39 85 d8 fb ff ff    	cmp    %eax,-0x428(%rbp)
   11483:	74 46                	je     114cb <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x7db>
        std::fprintf(stderr, "[%s:%d]: unmatched particle number %d vs %d\n", __FILE__, __LINE__, n, *nParticles);
   11485:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   1148c:	8b 08                	mov    (%rax),%ecx
   1148e:	48 8b 05 8b cd 00 00 	mov    0xcd8b(%rip),%rax        # 1e220 <stderr@GLIBC_2.2.5>
   11495:	8b 95 d8 fb ff ff    	mov    -0x428(%rbp),%edx
   1149b:	41 89 c9             	mov    %ecx,%r9d
   1149e:	41 89 d0             	mov    %edx,%r8d
   114a1:	b9 8b 00 00 00       	mov    $0x8b,%ecx
   114a6:	48 8d 15 f3 52 00 00 	lea    0x52f3(%rip),%rdx        # 167a0 <_ZStL19piecewise_construct+0x40>
   114ad:	48 8d 35 54 54 00 00 	lea    0x5454(%rip),%rsi        # 16908 <_ZStL19piecewise_construct+0x1a8>
   114b4:	48 89 c7             	mov    %rax,%rdi
   114b7:	b8 00 00 00 00       	mov    $0x0,%eax
   114bc:	e8 9f 92 ff ff       	call   a760 <fprintf@plt>
        exit(EXIT_FAILURE);
   114c1:	bf 01 00 00 00       	mov    $0x1,%edi
   114c6:	e8 75 90 ff ff       	call   a540 <exit@plt>
    }
}
   114cb:	48 8d 85 50 fc ff ff 	lea    -0x3b0(%rbp),%rax
   114d2:	48 89 c7             	mov    %rax,%rdi
   114d5:	e8 36 8d ff ff       	call   a210 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEED1Ev@plt>
   114da:	48 8d 85 30 fc ff ff 	lea    -0x3d0(%rbp),%rax
   114e1:	48 89 c7             	mov    %rax,%rdi
   114e4:	e8 37 8c ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   114e9:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   114f0:	48 89 c7             	mov    %rax,%rdi
   114f3:	e8 28 8a ff ff       	call   9f20 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEED1Ev@plt>
   114f8:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   114fc:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   11503:	00 00 
   11505:	74 70                	je     11577 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x887>
   11507:	eb 69                	jmp    11572 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x882>
   11509:	f3 0f 1e fa          	endbr64 
        std::fprintf(stderr, "[%s:%d]: can not open the file %s\n", __FILE__, __LINE__, filename);
   1150d:	48 89 c3             	mov    %rax,%rbx
   11510:	48 8d 85 50 fc ff ff 	lea    -0x3b0(%rbp),%rax
   11517:	48 89 c7             	mov    %rax,%rdi
   1151a:	e8 01 8c ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   1151f:	eb 37                	jmp    11558 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x868>
   11521:	f3 0f 1e fa          	endbr64 
}
   11525:	48 89 c3             	mov    %rax,%rbx
   11528:	48 8d 85 50 fc ff ff 	lea    -0x3b0(%rbp),%rax
   1152f:	48 89 c7             	mov    %rax,%rdi
   11532:	e8 d9 8c ff ff       	call   a210 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEED1Ev@plt>
   11537:	eb 07                	jmp    11540 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x850>
   11539:	f3 0f 1e fa          	endbr64 
   1153d:	48 89 c3             	mov    %rax,%rbx
   11540:	48 8d 85 30 fc ff ff 	lea    -0x3d0(%rbp),%rax
   11547:	48 89 c7             	mov    %rax,%rdi
   1154a:	e8 d1 8b ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   1154f:	eb 07                	jmp    11558 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x868>
   11551:	f3 0f 1e fa          	endbr64 
   11555:	48 89 c3             	mov    %rax,%rbx
   11558:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   1155f:	48 89 c7             	mov    %rax,%rdi
   11562:	e8 b9 89 ff ff       	call   9f20 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEED1Ev@plt>
   11567:	48 89 d8             	mov    %rbx,%rax
   1156a:	48 89 c7             	mov    %rax,%rdi
   1156d:	e8 3e 94 ff ff       	call   a9b0 <_Unwind_Resume@plt>
   11572:	e8 69 8f ff ff       	call   a4e0 <__stack_chk_fail@plt>
   11577:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   1157b:	c9                   	leave  
   1157c:	c3                   	ret    

000000000001157d <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti>:

void store_particles(const std::string filename, particle_t* particles, int nParticles) {
   1157d:	f3 0f 1e fa          	endbr64 
   11581:	55                   	push   %rbp
   11582:	48 89 e5             	mov    %rsp,%rbp
   11585:	53                   	push   %rbx
   11586:	48 81 ec 78 02 00 00 	sub    $0x278,%rsp
   1158d:	48 89 bd 98 fd ff ff 	mov    %rdi,-0x268(%rbp)
   11594:	48 89 b5 90 fd ff ff 	mov    %rsi,-0x270(%rbp)
   1159b:	89 95 8c fd ff ff    	mov    %edx,-0x274(%rbp)
   115a1:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   115a8:	00 00 
   115aa:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   115ae:	31 c0                	xor    %eax,%eax
    std::ofstream file(filename, std::ios::out|std::ios::trunc);
   115b0:	be 20 00 00 00       	mov    $0x20,%esi
   115b5:	bf 10 00 00 00       	mov    $0x10,%edi
   115ba:	e8 e7 11 00 00       	call   127a6 <_ZStorSt13_Ios_OpenmodeS_>
   115bf:	89 c2                	mov    %eax,%edx
   115c1:	48 8b 8d 98 fd ff ff 	mov    -0x268(%rbp),%rcx
   115c8:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   115cf:	48 89 ce             	mov    %rcx,%rsi
   115d2:	48 89 c7             	mov    %rax,%rdi
   115d5:	e8 16 95 ff ff       	call   aaf0 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEEC1ERKNSt7__cxx1112basic_stringIcS1_SaIcEEESt13_Ios_Openmode@plt>
    if(!file.is_open()) {
   115da:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   115e1:	48 89 c7             	mov    %rax,%rdi
   115e4:	e8 67 8c ff ff       	call   a250 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEE7is_openEv@plt>
   115e9:	83 f0 01             	xor    $0x1,%eax
   115ec:	84 c0                	test   %al,%al
   115ee:	74 63                	je     11653 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0xd6>
        std::fprintf(stderr, "[%s:%d]: can not open the file %s\n", __FILE__, __LINE__, filename);
   115f0:	48 8b 95 98 fd ff ff 	mov    -0x268(%rbp),%rdx
   115f7:	48 8d 85 c0 fd ff ff 	lea    -0x240(%rbp),%rax
   115fe:	48 89 d6             	mov    %rdx,%rsi
   11601:	48 89 c7             	mov    %rax,%rdi
   11604:	e8 d7 89 ff ff       	call   9fe0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
   11609:	48 8b 05 10 cc 00 00 	mov    0xcc10(%rip),%rax        # 1e220 <stderr@GLIBC_2.2.5>
   11610:	48 8d 95 c0 fd ff ff 	lea    -0x240(%rbp),%rdx
   11617:	49 89 d0             	mov    %rdx,%r8
   1161a:	b9 93 00 00 00       	mov    $0x93,%ecx
   1161f:	48 8d 15 7a 51 00 00 	lea    0x517a(%rip),%rdx        # 167a0 <_ZStL19piecewise_construct+0x40>
   11626:	48 8d 35 33 52 00 00 	lea    0x5233(%rip),%rsi        # 16860 <_ZStL19piecewise_construct+0x100>
   1162d:	48 89 c7             	mov    %rax,%rdi
   11630:	b8 00 00 00 00       	mov    $0x0,%eax
   11635:	e8 26 91 ff ff       	call   a760 <fprintf@plt>
   1163a:	48 8d 85 c0 fd ff ff 	lea    -0x240(%rbp),%rax
   11641:	48 89 c7             	mov    %rax,%rdi
   11644:	e8 d7 8a ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
        exit(EXIT_FAILURE);
   11649:	bf 01 00 00 00       	mov    $0x1,%edi
   1164e:	e8 ed 8e ff ff       	call   a540 <exit@plt>
    }
    int space_dimension = particles->ndim;
   11653:	48 8b 85 90 fd ff ff 	mov    -0x270(%rbp),%rax
   1165a:	8b 40 30             	mov    0x30(%rax),%eax
   1165d:	89 85 b4 fd ff ff    	mov    %eax,-0x24c(%rbp)
    file << nParticles << " " << space_dimension << std::endl;
   11663:	8b 95 8c fd ff ff    	mov    -0x274(%rbp),%edx
   11669:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11670:	89 d6                	mov    %edx,%esi
   11672:	48 89 c7             	mov    %rax,%rdi
   11675:	e8 06 93 ff ff       	call   a980 <_ZNSolsEi@plt>
   1167a:	48 89 c2             	mov    %rax,%rdx
   1167d:	48 8d 05 6e 52 00 00 	lea    0x526e(%rip),%rax        # 168f2 <_ZStL19piecewise_construct+0x192>
   11684:	48 89 c6             	mov    %rax,%rsi
   11687:	48 89 d7             	mov    %rdx,%rdi
   1168a:	e8 41 8d ff ff       	call   a3d0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
   1168f:	48 89 c2             	mov    %rax,%rdx
   11692:	8b 85 b4 fd ff ff    	mov    -0x24c(%rbp),%eax
   11698:	89 c6                	mov    %eax,%esi
   1169a:	48 89 d7             	mov    %rdx,%rdi
   1169d:	e8 de 92 ff ff       	call   a980 <_ZNSolsEi@plt>
   116a2:	48 8b 15 27 b9 00 00 	mov    0xb927(%rip),%rdx        # 1cfd0 <_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_@GLIBCXX_3.4>
   116a9:	48 89 d6             	mov    %rdx,%rsi
   116ac:	48 89 c7             	mov    %rax,%rdi
   116af:	e8 5c 8d ff ff       	call   a410 <_ZNSolsEPFRSoS_E@plt>
    file.flush();
   116b4:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   116bb:	48 89 c7             	mov    %rax,%rdi
   116be:	e8 cd 8a ff ff       	call   a190 <_ZNSo5flushEv@plt>
    particle_t* pp = particles;
   116c3:	48 8b 85 90 fd ff ff 	mov    -0x270(%rbp),%rax
   116ca:	48 89 85 b8 fd ff ff 	mov    %rax,-0x248(%rbp)
    for(int n=0; n<nParticles; n++) {
   116d1:	c7 85 a0 fd ff ff 00 	movl   $0x0,-0x260(%rbp)
   116d8:	00 00 00 
   116db:	e9 17 03 00 00       	jmp    119f7 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x47a>
        pp = particles + n;
   116e0:	8b 85 a0 fd ff ff    	mov    -0x260(%rbp),%eax
   116e6:	48 98                	cltq   
   116e8:	48 c1 e0 06          	shl    $0x6,%rax
   116ec:	48 89 c2             	mov    %rax,%rdx
   116ef:	48 8b 85 90 fd ff ff 	mov    -0x270(%rbp),%rax
   116f6:	48 01 d0             	add    %rdx,%rax
   116f9:	48 89 85 b8 fd ff ff 	mov    %rax,-0x248(%rbp)
        if (!pp->enabled) {
   11700:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11707:	0f b6 40 38          	movzbl 0x38(%rax),%eax
   1170b:	83 f0 01             	xor    $0x1,%eax
   1170e:	84 c0                	test   %al,%al
   11710:	74 2b                	je     1173d <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x1c0>
            std::cerr << "Before store, the data should already been shrinked and no disabled particle" << std::endl;
   11712:	48 8d 05 1f 52 00 00 	lea    0x521f(%rip),%rax        # 16938 <_ZStL19piecewise_construct+0x1d8>
   11719:	48 89 c6             	mov    %rax,%rsi
   1171c:	48 8d 05 1d cb 00 00 	lea    0xcb1d(%rip),%rax        # 1e240 <_ZSt4cerr@GLIBCXX_3.4>
   11723:	48 89 c7             	mov    %rax,%rdi
   11726:	e8 a5 8c ff ff       	call   a3d0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
   1172b:	48 8b 15 9e b8 00 00 	mov    0xb89e(%rip),%rdx        # 1cfd0 <_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_@GLIBCXX_3.4>
   11732:	48 89 d6             	mov    %rdx,%rsi
   11735:	48 89 c7             	mov    %rax,%rdi
   11738:	e8 d3 8c ff ff       	call   a410 <_ZNSolsEPFRSoS_E@plt>
            // assert(pp->enabled);
        }
        file << std::scientific << pp->id << " ";
   1173d:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11744:	48 8d 15 d7 10 00 00 	lea    0x10d7(%rip),%rdx        # 12822 <_ZSt10scientificRSt8ios_base>
   1174b:	48 89 d6             	mov    %rdx,%rsi
   1174e:	48 89 c7             	mov    %rax,%rdi
   11751:	e8 3a 88 ff ff       	call   9f90 <_ZNSolsEPFRSt8ios_baseS0_E@plt>
   11756:	48 89 c2             	mov    %rax,%rdx
   11759:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11760:	48 8b 00             	mov    (%rax),%rax
   11763:	48 89 c6             	mov    %rax,%rsi
   11766:	48 89 d7             	mov    %rdx,%rdi
   11769:	e8 82 87 ff ff       	call   9ef0 <_ZNSolsEm@plt>
   1176e:	48 89 c2             	mov    %rax,%rdx
   11771:	48 8d 05 7a 51 00 00 	lea    0x517a(%rip),%rax        # 168f2 <_ZStL19piecewise_construct+0x192>
   11778:	48 89 c6             	mov    %rax,%rsi
   1177b:	48 89 d7             	mov    %rdx,%rdi
   1177e:	e8 4d 8c ff ff       	call   a3d0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
        for(int nn=0; nn<(pp->ndim); nn++)
   11783:	c7 85 a4 fd ff ff 00 	movl   $0x0,-0x25c(%rbp)
   1178a:	00 00 00 
   1178d:	eb 62                	jmp    117f1 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x274>
            file << std::scientific << *((pp->position)+nn) << " ";
   1178f:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11796:	48 8d 15 85 10 00 00 	lea    0x1085(%rip),%rdx        # 12822 <_ZSt10scientificRSt8ios_base>
   1179d:	48 89 d6             	mov    %rdx,%rsi
   117a0:	48 89 c7             	mov    %rax,%rdi
   117a3:	e8 e8 87 ff ff       	call   9f90 <_ZNSolsEPFRSt8ios_baseS0_E@plt>
   117a8:	48 89 c2             	mov    %rax,%rdx
   117ab:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   117b2:	48 8b 48 08          	mov    0x8(%rax),%rcx
   117b6:	8b 85 a4 fd ff ff    	mov    -0x25c(%rbp),%eax
   117bc:	48 98                	cltq   
   117be:	48 c1 e0 03          	shl    $0x3,%rax
   117c2:	48 01 c8             	add    %rcx,%rax
   117c5:	48 8b 00             	mov    (%rax),%rax
   117c8:	66 48 0f 6e c0       	movq   %rax,%xmm0
   117cd:	48 89 d7             	mov    %rdx,%rdi
   117d0:	e8 ab 92 ff ff       	call   aa80 <_ZNSolsEd@plt>
   117d5:	48 89 c2             	mov    %rax,%rdx
   117d8:	48 8d 05 13 51 00 00 	lea    0x5113(%rip),%rax        # 168f2 <_ZStL19piecewise_construct+0x192>
   117df:	48 89 c6             	mov    %rax,%rsi
   117e2:	48 89 d7             	mov    %rdx,%rdi
   117e5:	e8 e6 8b ff ff       	call   a3d0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
        for(int nn=0; nn<(pp->ndim); nn++)
   117ea:	83 85 a4 fd ff ff 01 	addl   $0x1,-0x25c(%rbp)
   117f1:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   117f8:	8b 40 30             	mov    0x30(%rax),%eax
   117fb:	39 85 a4 fd ff ff    	cmp    %eax,-0x25c(%rbp)
   11801:	7c 8c                	jl     1178f <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x212>
        for(int nn=0; nn<(pp->ndim); nn++)
   11803:	c7 85 a8 fd ff ff 00 	movl   $0x0,-0x258(%rbp)
   1180a:	00 00 00 
   1180d:	eb 62                	jmp    11871 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x2f4>
            file << std::scientific << *((pp->velocity)+nn) << " ";
   1180f:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11816:	48 8d 15 05 10 00 00 	lea    0x1005(%rip),%rdx        # 12822 <_ZSt10scientificRSt8ios_base>
   1181d:	48 89 d6             	mov    %rdx,%rsi
   11820:	48 89 c7             	mov    %rax,%rdi
   11823:	e8 68 87 ff ff       	call   9f90 <_ZNSolsEPFRSt8ios_baseS0_E@plt>
   11828:	48 89 c2             	mov    %rax,%rdx
   1182b:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11832:	48 8b 48 10          	mov    0x10(%rax),%rcx
   11836:	8b 85 a8 fd ff ff    	mov    -0x258(%rbp),%eax
   1183c:	48 98                	cltq   
   1183e:	48 c1 e0 03          	shl    $0x3,%rax
   11842:	48 01 c8             	add    %rcx,%rax
   11845:	48 8b 00             	mov    (%rax),%rax
   11848:	66 48 0f 6e c0       	movq   %rax,%xmm0
   1184d:	48 89 d7             	mov    %rdx,%rdi
   11850:	e8 2b 92 ff ff       	call   aa80 <_ZNSolsEd@plt>
   11855:	48 89 c2             	mov    %rax,%rdx
   11858:	48 8d 05 93 50 00 00 	lea    0x5093(%rip),%rax        # 168f2 <_ZStL19piecewise_construct+0x192>
   1185f:	48 89 c6             	mov    %rax,%rsi
   11862:	48 89 d7             	mov    %rdx,%rdi
   11865:	e8 66 8b ff ff       	call   a3d0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
        for(int nn=0; nn<(pp->ndim); nn++)
   1186a:	83 85 a8 fd ff ff 01 	addl   $0x1,-0x258(%rbp)
   11871:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11878:	8b 40 30             	mov    0x30(%rax),%eax
   1187b:	39 85 a8 fd ff ff    	cmp    %eax,-0x258(%rbp)
   11881:	7c 8c                	jl     1180f <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x292>
        for(int nn=0; nn<(pp->ndim); nn++)
   11883:	c7 85 ac fd ff ff 00 	movl   $0x0,-0x254(%rbp)
   1188a:	00 00 00 
   1188d:	eb 62                	jmp    118f1 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x374>
            file << std::scientific << *((pp->acceleration)+nn) << " ";
   1188f:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11896:	48 8d 15 85 0f 00 00 	lea    0xf85(%rip),%rdx        # 12822 <_ZSt10scientificRSt8ios_base>
   1189d:	48 89 d6             	mov    %rdx,%rsi
   118a0:	48 89 c7             	mov    %rax,%rdi
   118a3:	e8 e8 86 ff ff       	call   9f90 <_ZNSolsEPFRSt8ios_baseS0_E@plt>
   118a8:	48 89 c2             	mov    %rax,%rdx
   118ab:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   118b2:	48 8b 48 18          	mov    0x18(%rax),%rcx
   118b6:	8b 85 ac fd ff ff    	mov    -0x254(%rbp),%eax
   118bc:	48 98                	cltq   
   118be:	48 c1 e0 03          	shl    $0x3,%rax
   118c2:	48 01 c8             	add    %rcx,%rax
   118c5:	48 8b 00             	mov    (%rax),%rax
   118c8:	66 48 0f 6e c0       	movq   %rax,%xmm0
   118cd:	48 89 d7             	mov    %rdx,%rdi
   118d0:	e8 ab 91 ff ff       	call   aa80 <_ZNSolsEd@plt>
   118d5:	48 89 c2             	mov    %rax,%rdx
   118d8:	48 8d 05 13 50 00 00 	lea    0x5013(%rip),%rax        # 168f2 <_ZStL19piecewise_construct+0x192>
   118df:	48 89 c6             	mov    %rax,%rsi
   118e2:	48 89 d7             	mov    %rdx,%rdi
   118e5:	e8 e6 8a ff ff       	call   a3d0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
        for(int nn=0; nn<(pp->ndim); nn++)
   118ea:	83 85 ac fd ff ff 01 	addl   $0x1,-0x254(%rbp)
   118f1:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   118f8:	8b 40 30             	mov    0x30(%rax),%eax
   118fb:	39 85 ac fd ff ff    	cmp    %eax,-0x254(%rbp)
   11901:	7c 8c                	jl     1188f <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x312>
        file << pp->nfeat << " ";
   11903:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   1190a:	8b 50 34             	mov    0x34(%rax),%edx
   1190d:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11914:	89 d6                	mov    %edx,%esi
   11916:	48 89 c7             	mov    %rax,%rdi
   11919:	e8 62 90 ff ff       	call   a980 <_ZNSolsEi@plt>
   1191e:	48 89 c2             	mov    %rax,%rdx
   11921:	48 8d 05 ca 4f 00 00 	lea    0x4fca(%rip),%rax        # 168f2 <_ZStL19piecewise_construct+0x192>
   11928:	48 89 c6             	mov    %rax,%rsi
   1192b:	48 89 d7             	mov    %rdx,%rdi
   1192e:	e8 9d 8a ff ff       	call   a3d0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
        for(int nn=0; nn<pp->nfeat; nn++) {
   11933:	c7 85 b0 fd ff ff 00 	movl   $0x0,-0x250(%rbp)
   1193a:	00 00 00 
   1193d:	e9 98 00 00 00       	jmp    119da <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x45d>
            file << std::scientific << *(pp->features+nn);
   11942:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11949:	48 8d 15 d2 0e 00 00 	lea    0xed2(%rip),%rdx        # 12822 <_ZSt10scientificRSt8ios_base>
   11950:	48 89 d6             	mov    %rdx,%rsi
   11953:	48 89 c7             	mov    %rax,%rdi
   11956:	e8 35 86 ff ff       	call   9f90 <_ZNSolsEPFRSt8ios_baseS0_E@plt>
   1195b:	48 89 c2             	mov    %rax,%rdx
   1195e:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11965:	48 8b 48 20          	mov    0x20(%rax),%rcx
   11969:	8b 85 b0 fd ff ff    	mov    -0x250(%rbp),%eax
   1196f:	48 98                	cltq   
   11971:	48 c1 e0 03          	shl    $0x3,%rax
   11975:	48 01 c8             	add    %rcx,%rax
   11978:	48 8b 00             	mov    (%rax),%rax
   1197b:	66 48 0f 6e c0       	movq   %rax,%xmm0
   11980:	48 89 d7             	mov    %rdx,%rdi
   11983:	e8 f8 90 ff ff       	call   aa80 <_ZNSolsEd@plt>
            if(nn+1 != pp->nfeat)
   11988:	8b 85 b0 fd ff ff    	mov    -0x250(%rbp),%eax
   1198e:	8d 50 01             	lea    0x1(%rax),%edx
   11991:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11998:	8b 40 34             	mov    0x34(%rax),%eax
   1199b:	39 c2                	cmp    %eax,%edx
   1199d:	74 1b                	je     119ba <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x43d>
                file << " ";
   1199f:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   119a6:	48 8d 15 45 4f 00 00 	lea    0x4f45(%rip),%rdx        # 168f2 <_ZStL19piecewise_construct+0x192>
   119ad:	48 89 d6             	mov    %rdx,%rsi
   119b0:	48 89 c7             	mov    %rax,%rdi
   119b3:	e8 18 8a ff ff       	call   a3d0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
   119b8:	eb 19                	jmp    119d3 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x456>
            else
                file << std::endl;
   119ba:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   119c1:	48 8b 15 08 b6 00 00 	mov    0xb608(%rip),%rdx        # 1cfd0 <_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_@GLIBCXX_3.4>
   119c8:	48 89 d6             	mov    %rdx,%rsi
   119cb:	48 89 c7             	mov    %rax,%rdi
   119ce:	e8 3d 8a ff ff       	call   a410 <_ZNSolsEPFRSoS_E@plt>
        for(int nn=0; nn<pp->nfeat; nn++) {
   119d3:	83 85 b0 fd ff ff 01 	addl   $0x1,-0x250(%rbp)
   119da:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   119e1:	8b 40 34             	mov    0x34(%rax),%eax
   119e4:	39 85 b0 fd ff ff    	cmp    %eax,-0x250(%rbp)
   119ea:	0f 8c 52 ff ff ff    	jl     11942 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x3c5>
    for(int n=0; n<nParticles; n++) {
   119f0:	83 85 a0 fd ff ff 01 	addl   $0x1,-0x260(%rbp)
   119f7:	8b 85 a0 fd ff ff    	mov    -0x260(%rbp),%eax
   119fd:	3b 85 8c fd ff ff    	cmp    -0x274(%rbp),%eax
   11a03:	0f 8c d7 fc ff ff    	jl     116e0 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x163>
        }
    }
    file.close();
   11a09:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11a10:	48 89 c7             	mov    %rax,%rdi
   11a13:	e8 08 8d ff ff       	call   a720 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEE5closeEv@plt>
}
   11a18:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11a1f:	48 89 c7             	mov    %rax,%rdi
   11a22:	e8 59 8d ff ff       	call   a780 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEED1Ev@plt>
   11a27:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   11a2b:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   11a32:	00 00 
   11a34:	74 40                	je     11a76 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x4f9>
   11a36:	eb 39                	jmp    11a71 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x4f4>
   11a38:	f3 0f 1e fa          	endbr64 
        std::fprintf(stderr, "[%s:%d]: can not open the file %s\n", __FILE__, __LINE__, filename);
   11a3c:	48 89 c3             	mov    %rax,%rbx
   11a3f:	48 8d 85 c0 fd ff ff 	lea    -0x240(%rbp),%rax
   11a46:	48 89 c7             	mov    %rax,%rdi
   11a49:	e8 d2 86 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   11a4e:	eb 07                	jmp    11a57 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x4da>
   11a50:	f3 0f 1e fa          	endbr64 
}
   11a54:	48 89 c3             	mov    %rax,%rbx
   11a57:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11a5e:	48 89 c7             	mov    %rax,%rdi
   11a61:	e8 1a 8d ff ff       	call   a780 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEED1Ev@plt>
   11a66:	48 89 d8             	mov    %rbx,%rax
   11a69:	48 89 c7             	mov    %rax,%rdi
   11a6c:	e8 3f 8f ff ff       	call   a9b0 <_Unwind_Resume@plt>
   11a71:	e8 6a 8a ff ff       	call   a4e0 <__stack_chk_fail@plt>
   11a76:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   11a7a:	c9                   	leave  
   11a7b:	c3                   	ret    

0000000000011a7c <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t>:

void load_particles(const std::string filename, chunk_particles_t** particleChunk) {
   11a7c:	f3 0f 1e fa          	endbr64 
   11a80:	55                   	push   %rbp
   11a81:	48 89 e5             	mov    %rsp,%rbp
   11a84:	53                   	push   %rbx
   11a85:	48 83 ec 58          	sub    $0x58,%rsp
   11a89:	48 89 7d a8          	mov    %rdi,-0x58(%rbp)
   11a8d:	48 89 75 a0          	mov    %rsi,-0x60(%rbp)
   11a91:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   11a98:	00 00 
   11a9a:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   11a9e:	31 c0                	xor    %eax,%eax
    if(*particleChunk == nullptr)
   11aa0:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   11aa4:	48 8b 00             	mov    (%rax),%rax
   11aa7:	48 85 c0             	test   %rax,%rax
   11aaa:	75 16                	jne    11ac2 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t+0x46>
        *particleChunk = static_cast<chunk_particles_t*>(malloc(sizeof(chunk_particles_t)));
   11aac:	bf 10 00 00 00       	mov    $0x10,%edi
   11ab1:	e8 7a 8b ff ff       	call   a630 <malloc@plt>
   11ab6:	48 89 c2             	mov    %rax,%rdx
   11ab9:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   11abd:	48 89 10             	mov    %rdx,(%rax)
   11ac0:	eb 36                	jmp    11af8 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t+0x7c>
    else
        assert((*particleChunk)->nParticle == 0);
   11ac2:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   11ac6:	48 8b 00             	mov    (%rax),%rax
   11ac9:	8b 40 08             	mov    0x8(%rax),%eax
   11acc:	85 c0                	test   %eax,%eax
   11ace:	74 28                	je     11af8 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t+0x7c>
   11ad0:	48 8d 05 b1 4e 00 00 	lea    0x4eb1(%rip),%rax        # 16988 <_ZStL19piecewise_construct+0x228>
   11ad7:	48 89 c1             	mov    %rax,%rcx
   11ada:	ba b7 00 00 00       	mov    $0xb7,%edx
   11adf:	48 8d 05 ba 4c 00 00 	lea    0x4cba(%rip),%rax        # 167a0 <_ZStL19piecewise_construct+0x40>
   11ae6:	48 89 c6             	mov    %rax,%rsi
   11ae9:	48 8d 05 d0 4e 00 00 	lea    0x4ed0(%rip),%rax        # 169c0 <_ZStL19piecewise_construct+0x260>
   11af0:	48 89 c7             	mov    %rax,%rdi
   11af3:	e8 58 86 ff ff       	call   a150 <__assert_fail@plt>
    particle_t* particles;
    int nParticles;
    load_particles(filename, &particles, &nParticles);
   11af8:	48 8b 55 a8          	mov    -0x58(%rbp),%rdx
   11afc:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11b00:	48 89 d6             	mov    %rdx,%rsi
   11b03:	48 89 c7             	mov    %rax,%rdi
   11b06:	e8 d5 84 ff ff       	call   9fe0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
   11b0b:	48 8d 55 b4          	lea    -0x4c(%rbp),%rdx
   11b0f:	48 8d 4d b8          	lea    -0x48(%rbp),%rcx
   11b13:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11b17:	48 89 ce             	mov    %rcx,%rsi
   11b1a:	48 89 c7             	mov    %rax,%rdi
   11b1d:	e8 ce f1 ff ff       	call   10cf0 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi>
   11b22:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11b26:	48 89 c7             	mov    %rax,%rdi
   11b29:	e8 f2 85 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    (*particleChunk)->particles = particles;
   11b2e:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   11b32:	48 8b 00             	mov    (%rax),%rax
   11b35:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
   11b39:	48 89 10             	mov    %rdx,(%rax)
    (*particleChunk)->nParticle = nParticles;
   11b3c:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   11b40:	48 8b 00             	mov    (%rax),%rax
   11b43:	8b 55 b4             	mov    -0x4c(%rbp),%edx
   11b46:	89 50 08             	mov    %edx,0x8(%rax)
}
   11b49:	90                   	nop
   11b4a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   11b4e:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   11b55:	00 00 
   11b57:	74 25                	je     11b7e <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t+0x102>
   11b59:	eb 1e                	jmp    11b79 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t+0xfd>
   11b5b:	f3 0f 1e fa          	endbr64 
    load_particles(filename, &particles, &nParticles);
   11b5f:	48 89 c3             	mov    %rax,%rbx
   11b62:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11b66:	48 89 c7             	mov    %rax,%rdi
   11b69:	e8 b2 85 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   11b6e:	48 89 d8             	mov    %rbx,%rax
   11b71:	48 89 c7             	mov    %rax,%rdi
   11b74:	e8 37 8e ff ff       	call   a9b0 <_Unwind_Resume@plt>
}
   11b79:	e8 62 89 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   11b7e:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   11b82:	c9                   	leave  
   11b83:	c3                   	ret    

0000000000011b84 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t>:

void store_particles(const std::string filename, chunk_particles_t* particleChunk) {
   11b84:	f3 0f 1e fa          	endbr64 
   11b88:	55                   	push   %rbp
   11b89:	48 89 e5             	mov    %rsp,%rbp
   11b8c:	41 54                	push   %r12
   11b8e:	53                   	push   %rbx
   11b8f:	48 83 ec 40          	sub    $0x40,%rsp
   11b93:	48 89 7d b8          	mov    %rdi,-0x48(%rbp)
   11b97:	48 89 75 b0          	mov    %rsi,-0x50(%rbp)
   11b9b:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   11ba2:	00 00 
   11ba4:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   11ba8:	31 c0                	xor    %eax,%eax
    store_particles(filename, particleChunk->particles, particleChunk->nParticle);
   11baa:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   11bae:	44 8b 60 08          	mov    0x8(%rax),%r12d
   11bb2:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   11bb6:	48 8b 18             	mov    (%rax),%rbx
   11bb9:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
   11bbd:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11bc1:	48 89 d6             	mov    %rdx,%rsi
   11bc4:	48 89 c7             	mov    %rax,%rdi
   11bc7:	e8 14 84 ff ff       	call   9fe0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
   11bcc:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11bd0:	44 89 e2             	mov    %r12d,%edx
   11bd3:	48 89 de             	mov    %rbx,%rsi
   11bd6:	48 89 c7             	mov    %rax,%rdi
   11bd9:	e8 9f f9 ff ff       	call   1157d <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti>
   11bde:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11be2:	48 89 c7             	mov    %rax,%rdi
   11be5:	e8 36 85 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
}
   11bea:	90                   	nop
   11beb:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   11bef:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   11bf6:	00 00 
   11bf8:	74 25                	je     11c1f <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t+0x9b>
   11bfa:	eb 1e                	jmp    11c1a <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t+0x96>
   11bfc:	f3 0f 1e fa          	endbr64 
    store_particles(filename, particleChunk->particles, particleChunk->nParticle);
   11c00:	48 89 c3             	mov    %rax,%rbx
   11c03:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11c07:	48 89 c7             	mov    %rax,%rdi
   11c0a:	e8 11 85 ff ff       	call   a120 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   11c0f:	48 89 d8             	mov    %rbx,%rax
   11c12:	48 89 c7             	mov    %rax,%rdi
   11c15:	e8 96 8d ff ff       	call   a9b0 <_Unwind_Resume@plt>
}
   11c1a:	e8 c1 88 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   11c1f:	48 83 c4 40          	add    $0x40,%rsp
   11c23:	5b                   	pop    %rbx
   11c24:	41 5c                	pop    %r12
   11c26:	5d                   	pop    %rbp
   11c27:	c3                   	ret    

0000000000011c28 <_Z14free_particlesP10particle_t>:

void free_particles(particle_t* particles) {
   11c28:	f3 0f 1e fa          	endbr64 
   11c2c:	55                   	push   %rbp
   11c2d:	48 89 e5             	mov    %rsp,%rbp
   11c30:	48 83 ec 10          	sub    $0x10,%rsp
   11c34:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    free(particles->position);
   11c38:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11c3c:	48 8b 40 08          	mov    0x8(%rax),%rax
   11c40:	48 89 c7             	mov    %rax,%rdi
   11c43:	e8 b8 88 ff ff       	call   a500 <free@plt>
    free(particles->velocity);
   11c48:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11c4c:	48 8b 40 10          	mov    0x10(%rax),%rax
   11c50:	48 89 c7             	mov    %rax,%rdi
   11c53:	e8 a8 88 ff ff       	call   a500 <free@plt>
    free(particles->acceleration);
   11c58:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11c5c:	48 8b 40 18          	mov    0x18(%rax),%rax
   11c60:	48 89 c7             	mov    %rax,%rdi
   11c63:	e8 98 88 ff ff       	call   a500 <free@plt>
    free(particles->features);
   11c68:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11c6c:	48 8b 40 20          	mov    0x20(%rax),%rax
   11c70:	48 89 c7             	mov    %rax,%rdi
   11c73:	e8 88 88 ff ff       	call   a500 <free@plt>
    free(particles);
   11c78:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11c7c:	48 89 c7             	mov    %rax,%rdi
   11c7f:	e8 7c 88 ff ff       	call   a500 <free@plt>
}
   11c84:	90                   	nop
   11c85:	c9                   	leave  
   11c86:	c3                   	ret    

0000000000011c87 <_Z14free_particlesP17chunk_particles_t>:

void free_particles(chunk_particles_t* particles) {
   11c87:	f3 0f 1e fa          	endbr64 
   11c8b:	55                   	push   %rbp
   11c8c:	48 89 e5             	mov    %rsp,%rbp
   11c8f:	48 83 ec 10          	sub    $0x10,%rsp
   11c93:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    free_particles(particles->particles);
   11c97:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11c9b:	48 8b 00             	mov    (%rax),%rax
   11c9e:	48 89 c7             	mov    %rax,%rdi
   11ca1:	e8 82 ff ff ff       	call   11c28 <_Z14free_particlesP10particle_t>
    particles->particles = nullptr;
   11ca6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11caa:	48 c7 00 00 00 00 00 	movq   $0x0,(%rax)
    particles->nParticle = 0;
   11cb1:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11cb5:	c7 40 08 00 00 00 00 	movl   $0x0,0x8(%rax)
}
   11cbc:	90                   	nop
   11cbd:	c9                   	leave  
   11cbe:	c3                   	ret    

0000000000011cbf <_Z18clearAccelerationsP10particle_ti>:

void clearAccelerations(particle_t* particles, int nParticle) {
   11cbf:	f3 0f 1e fa          	endbr64 
   11cc3:	55                   	push   %rbp
   11cc4:	48 89 e5             	mov    %rsp,%rbp
   11cc7:	48 83 ec 40          	sub    $0x40,%rsp
   11ccb:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
   11ccf:	89 75 c4             	mov    %esi,-0x3c(%rbp)
   11cd2:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   11cd9:	00 00 
   11cdb:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   11cdf:	31 c0                	xor    %eax,%eax
    particle_t* pt = particles;
   11ce1:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   11ce5:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
    #pragma omp parallel for
   11ce9:	8b 45 c4             	mov    -0x3c(%rbp),%eax
   11cec:	89 45 e8             	mov    %eax,-0x18(%rbp)
   11cef:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   11cf3:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   11cf7:	48 8d 45 e0          	lea    -0x20(%rbp),%rax
   11cfb:	b9 00 00 00 00       	mov    $0x0,%ecx
   11d00:	ba 00 00 00 00       	mov    $0x0,%edx
   11d05:	48 89 c6             	mov    %rax,%rsi
   11d08:	48 8d 05 5f 04 00 00 	lea    0x45f(%rip),%rax        # 1216e <_Z18clearAccelerationsP10particle_ti._omp_fn.0>
   11d0f:	48 89 c7             	mov    %rax,%rdi
   11d12:	e8 99 8a ff ff       	call   a7b0 <GOMP_parallel@plt>
   11d17:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   11d1b:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
        #pragma UNROLL
        for(int dim=0; dim<pt->ndim; dim++)
            *(acceleration++) = 0;
        pt++;
    }
}
   11d1f:	90                   	nop
   11d20:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11d24:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   11d2b:	00 00 
   11d2d:	74 05                	je     11d34 <_Z18clearAccelerationsP10particle_ti+0x75>
   11d2f:	e8 ac 87 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   11d34:	c9                   	leave  
   11d35:	c3                   	ret    

0000000000011d36 <_Z18clearAccelerationsP17chunk_particles_t>:

void clearAccelerations(chunk_particles_t* particlesChunk) {
   11d36:	f3 0f 1e fa          	endbr64 
   11d3a:	55                   	push   %rbp
   11d3b:	48 89 e5             	mov    %rsp,%rbp
   11d3e:	48 83 ec 10          	sub    $0x10,%rsp
   11d42:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    clearAccelerations(particlesChunk->particles, particlesChunk->nParticle);
   11d46:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11d4a:	8b 50 08             	mov    0x8(%rax),%edx
   11d4d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11d51:	48 8b 00             	mov    (%rax),%rax
   11d54:	89 d6                	mov    %edx,%esi
   11d56:	48 89 c7             	mov    %rax,%rdi
   11d59:	e8 61 ff ff ff       	call   11cbf <_Z18clearAccelerationsP10particle_ti>
}
   11d5e:	90                   	nop
   11d5f:	c9                   	leave  
   11d60:	c3                   	ret    

0000000000011d61 <_Z16update_particlesP10particle_tid>:

void update_particles(particle_t* particles, int nParticle, double timeStep) {
   11d61:	f3 0f 1e fa          	endbr64 
   11d65:	55                   	push   %rbp
   11d66:	48 89 e5             	mov    %rsp,%rbp
   11d69:	48 83 ec 40          	sub    $0x40,%rsp
   11d6d:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   11d71:	89 75 d4             	mov    %esi,-0x2c(%rbp)
   11d74:	f2 0f 11 45 c8       	movsd  %xmm0,-0x38(%rbp)
   11d79:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   11d80:	00 00 
   11d82:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   11d86:	31 c0                	xor    %eax,%eax
    clearAccelerations(particles, nParticle);
   11d88:	8b 55 d4             	mov    -0x2c(%rbp),%edx
   11d8b:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   11d8f:	89 d6                	mov    %edx,%esi
   11d91:	48 89 c7             	mov    %rax,%rdi
   11d94:	e8 26 ff ff ff       	call   11cbf <_Z18clearAccelerationsP10particle_ti>
    #pragma omp parallel for collapse(2)
   11d99:	8b 45 d4             	mov    -0x2c(%rbp),%eax
   11d9c:	89 45 e8             	mov    %eax,-0x18(%rbp)
   11d9f:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   11da3:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   11da7:	48 8d 45 e0          	lea    -0x20(%rbp),%rax
   11dab:	b9 00 00 00 00       	mov    $0x0,%ecx
   11db0:	ba 00 00 00 00       	mov    $0x0,%edx
   11db5:	48 89 c6             	mov    %rax,%rsi
   11db8:	48 8d 05 74 04 00 00 	lea    0x474(%rip),%rax        # 12233 <_Z16update_particlesP10particle_tid._omp_fn.0>
   11dbf:	48 89 c7             	mov    %rax,%rdi
   11dc2:	e8 e9 89 ff ff       	call   a7b0 <GOMP_parallel@plt>
    for(int i=0; i<nParticle; i++) 
        for(int j=0; j<i; j++) 
            updateAcceleration(particles+i, particles+j);

    #pragma omp parallel for
   11dc7:	f2 0f 10 45 c8       	movsd  -0x38(%rbp),%xmm0
   11dcc:	f2 0f 11 45 e8       	movsd  %xmm0,-0x18(%rbp)
   11dd1:	8b 45 d4             	mov    -0x2c(%rbp),%eax
   11dd4:	89 45 f0             	mov    %eax,-0x10(%rbp)
   11dd7:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   11ddb:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   11ddf:	48 8d 45 e0          	lea    -0x20(%rbp),%rax
   11de3:	b9 00 00 00 00       	mov    $0x0,%ecx
   11de8:	ba 00 00 00 00       	mov    $0x0,%edx
   11ded:	48 89 c6             	mov    %rax,%rsi
   11df0:	48 8d 05 a0 06 00 00 	lea    0x6a0(%rip),%rax        # 12497 <_Z16update_particlesP10particle_tid._omp_fn.1>
   11df7:	48 89 c7             	mov    %rax,%rdi
   11dfa:	e8 b1 89 ff ff       	call   a7b0 <GOMP_parallel@plt>
    for(int i=0; i<nParticle; i++) {
        updateVelocity(particles+i, timeStep);
        updatePosition(particles+i, timeStep);
    }
}
   11dff:	90                   	nop
   11e00:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11e04:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   11e0b:	00 00 
   11e0d:	74 05                	je     11e14 <_Z16update_particlesP10particle_tid+0xb3>
   11e0f:	e8 cc 86 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   11e14:	c9                   	leave  
   11e15:	c3                   	ret    

0000000000011e16 <_Z16update_particlesP17chunk_particles_td>:

void update_particles(chunk_particles_t* particlesChunk, double timeStep) {
   11e16:	f3 0f 1e fa          	endbr64 
   11e1a:	55                   	push   %rbp
   11e1b:	48 89 e5             	mov    %rsp,%rbp
   11e1e:	48 83 ec 10          	sub    $0x10,%rsp
   11e22:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   11e26:	f2 0f 11 45 f0       	movsd  %xmm0,-0x10(%rbp)
    update_particles(particlesChunk->particles, particlesChunk->nParticle, timeStep);
   11e2b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11e2f:	8b 50 08             	mov    0x8(%rax),%edx
   11e32:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11e36:	48 8b 00             	mov    (%rax),%rax
   11e39:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
   11e3d:	66 48 0f 6e c1       	movq   %rcx,%xmm0
   11e42:	89 d6                	mov    %edx,%esi
   11e44:	48 89 c7             	mov    %rax,%rdi
   11e47:	e8 15 ff ff ff       	call   11d61 <_Z16update_particlesP10particle_tid>
}
   11e4c:	90                   	nop
   11e4d:	c9                   	leave  
   11e4e:	c3                   	ret    

0000000000011e4f <_Z15alloc_particlesPP10particle_tiii>:

void alloc_particles(particle_t** particles, int nParticle, int ndim, int nfeat) {
   11e4f:	f3 0f 1e fa          	endbr64 
   11e53:	55                   	push   %rbp
   11e54:	48 89 e5             	mov    %rsp,%rbp
   11e57:	48 83 c4 80          	add    $0xffffffffffffff80,%rsp
   11e5b:	48 89 7d 98          	mov    %rdi,-0x68(%rbp)
   11e5f:	89 75 94             	mov    %esi,-0x6c(%rbp)
   11e62:	89 55 90             	mov    %edx,-0x70(%rbp)
   11e65:	89 4d 8c             	mov    %ecx,-0x74(%rbp)
   11e68:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   11e6f:	00 00 
   11e71:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   11e75:	31 c0                	xor    %eax,%eax
    assert(*particles == nullptr);
   11e77:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   11e7b:	48 8b 00             	mov    (%rax),%rax
   11e7e:	48 85 c0             	test   %rax,%rax
   11e81:	74 28                	je     11eab <_Z15alloc_particlesPP10particle_tiii+0x5c>
   11e83:	48 8d 05 5e 4b 00 00 	lea    0x4b5e(%rip),%rax        # 169e8 <_ZStL19piecewise_construct+0x288>
   11e8a:	48 89 c1             	mov    %rax,%rcx
   11e8d:	ba f4 00 00 00       	mov    $0xf4,%edx
   11e92:	48 8d 05 07 49 00 00 	lea    0x4907(%rip),%rax        # 167a0 <_ZStL19piecewise_construct+0x40>
   11e99:	48 89 c6             	mov    %rax,%rsi
   11e9c:	48 8d 05 77 4b 00 00 	lea    0x4b77(%rip),%rax        # 16a1a <_ZStL19piecewise_construct+0x2ba>
   11ea3:	48 89 c7             	mov    %rax,%rdi
   11ea6:	e8 a5 82 ff ff       	call   a150 <__assert_fail@plt>
    *particles = static_cast<particle_t*>(malloc(nParticle*sizeof(particle_t)));
   11eab:	8b 45 94             	mov    -0x6c(%rbp),%eax
   11eae:	48 98                	cltq   
   11eb0:	48 c1 e0 06          	shl    $0x6,%rax
   11eb4:	48 89 c7             	mov    %rax,%rdi
   11eb7:	e8 74 87 ff ff       	call   a630 <malloc@plt>
   11ebc:	48 89 c2             	mov    %rax,%rdx
   11ebf:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   11ec3:	48 89 10             	mov    %rdx,(%rax)
    double* pos = static_cast<double*>(malloc(nParticle*ndim*sizeof(double)));
   11ec6:	8b 45 94             	mov    -0x6c(%rbp),%eax
   11ec9:	0f af 45 90          	imul   -0x70(%rbp),%eax
   11ecd:	48 98                	cltq   
   11ecf:	48 c1 e0 03          	shl    $0x3,%rax
   11ed3:	48 89 c7             	mov    %rax,%rdi
   11ed6:	e8 55 87 ff ff       	call   a630 <malloc@plt>
   11edb:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
    double* vel = static_cast<double*>(malloc(nParticle*ndim*sizeof(double)));
   11edf:	8b 45 94             	mov    -0x6c(%rbp),%eax
   11ee2:	0f af 45 90          	imul   -0x70(%rbp),%eax
   11ee6:	48 98                	cltq   
   11ee8:	48 c1 e0 03          	shl    $0x3,%rax
   11eec:	48 89 c7             	mov    %rax,%rdi
   11eef:	e8 3c 87 ff ff       	call   a630 <malloc@plt>
   11ef4:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
    double* acc = static_cast<double*>(malloc(nParticle*ndim*sizeof(double)));
   11ef8:	8b 45 94             	mov    -0x6c(%rbp),%eax
   11efb:	0f af 45 90          	imul   -0x70(%rbp),%eax
   11eff:	48 98                	cltq   
   11f01:	48 c1 e0 03          	shl    $0x3,%rax
   11f05:	48 89 c7             	mov    %rax,%rdi
   11f08:	e8 23 87 ff ff       	call   a630 <malloc@plt>
   11f0d:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
    double* feat = static_cast<double*>(malloc(nParticle*nfeat*sizeof(double)));
   11f11:	8b 45 94             	mov    -0x6c(%rbp),%eax
   11f14:	0f af 45 8c          	imul   -0x74(%rbp),%eax
   11f18:	48 98                	cltq   
   11f1a:	48 c1 e0 03          	shl    $0x3,%rax
   11f1e:	48 89 c7             	mov    %rax,%rdi
   11f21:	e8 0a 87 ff ff       	call   a630 <malloc@plt>
   11f26:	48 89 45 b8          	mov    %rax,-0x48(%rbp)
    assert(*particles != nullptr);
   11f2a:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   11f2e:	48 8b 00             	mov    (%rax),%rax
   11f31:	48 85 c0             	test   %rax,%rax
   11f34:	75 28                	jne    11f5e <_Z15alloc_particlesPP10particle_tiii+0x10f>
   11f36:	48 8d 05 ab 4a 00 00 	lea    0x4aab(%rip),%rax        # 169e8 <_ZStL19piecewise_construct+0x288>
   11f3d:	48 89 c1             	mov    %rax,%rcx
   11f40:	ba fa 00 00 00       	mov    $0xfa,%edx
   11f45:	48 8d 05 54 48 00 00 	lea    0x4854(%rip),%rax        # 167a0 <_ZStL19piecewise_construct+0x40>
   11f4c:	48 89 c6             	mov    %rax,%rsi
   11f4f:	48 8d 05 da 4a 00 00 	lea    0x4ada(%rip),%rax        # 16a30 <_ZStL19piecewise_construct+0x2d0>
   11f56:	48 89 c7             	mov    %rax,%rdi
   11f59:	e8 f2 81 ff ff       	call   a150 <__assert_fail@plt>
    assert(pos != nullptr);     assert(vel != nullptr);
   11f5e:	48 83 7d a0 00       	cmpq   $0x0,-0x60(%rbp)
   11f63:	75 28                	jne    11f8d <_Z15alloc_particlesPP10particle_tiii+0x13e>
   11f65:	48 8d 05 7c 4a 00 00 	lea    0x4a7c(%rip),%rax        # 169e8 <_ZStL19piecewise_construct+0x288>
   11f6c:	48 89 c1             	mov    %rax,%rcx
   11f6f:	ba fb 00 00 00       	mov    $0xfb,%edx
   11f74:	48 8d 05 25 48 00 00 	lea    0x4825(%rip),%rax        # 167a0 <_ZStL19piecewise_construct+0x40>
   11f7b:	48 89 c6             	mov    %rax,%rsi
   11f7e:	48 8d 05 40 49 00 00 	lea    0x4940(%rip),%rax        # 168c5 <_ZStL19piecewise_construct+0x165>
   11f85:	48 89 c7             	mov    %rax,%rdi
   11f88:	e8 c3 81 ff ff       	call   a150 <__assert_fail@plt>
   11f8d:	48 83 7d a8 00       	cmpq   $0x0,-0x58(%rbp)
   11f92:	75 28                	jne    11fbc <_Z15alloc_particlesPP10particle_tiii+0x16d>
   11f94:	48 8d 05 4d 4a 00 00 	lea    0x4a4d(%rip),%rax        # 169e8 <_ZStL19piecewise_construct+0x288>
   11f9b:	48 89 c1             	mov    %rax,%rcx
   11f9e:	ba fb 00 00 00       	mov    $0xfb,%edx
   11fa3:	48 8d 05 f6 47 00 00 	lea    0x47f6(%rip),%rax        # 167a0 <_ZStL19piecewise_construct+0x40>
   11faa:	48 89 c6             	mov    %rax,%rsi
   11fad:	48 8d 05 20 49 00 00 	lea    0x4920(%rip),%rax        # 168d4 <_ZStL19piecewise_construct+0x174>
   11fb4:	48 89 c7             	mov    %rax,%rdi
   11fb7:	e8 94 81 ff ff       	call   a150 <__assert_fail@plt>
    assert(acc != nullptr);     assert(feat != nullptr);
   11fbc:	48 83 7d b0 00       	cmpq   $0x0,-0x50(%rbp)
   11fc1:	75 28                	jne    11feb <_Z15alloc_particlesPP10particle_tiii+0x19c>
   11fc3:	48 8d 05 1e 4a 00 00 	lea    0x4a1e(%rip),%rax        # 169e8 <_ZStL19piecewise_construct+0x288>
   11fca:	48 89 c1             	mov    %rax,%rcx
   11fcd:	ba fc 00 00 00       	mov    $0xfc,%edx
   11fd2:	48 8d 05 c7 47 00 00 	lea    0x47c7(%rip),%rax        # 167a0 <_ZStL19piecewise_construct+0x40>
   11fd9:	48 89 c6             	mov    %rax,%rsi
   11fdc:	48 8d 05 00 49 00 00 	lea    0x4900(%rip),%rax        # 168e3 <_ZStL19piecewise_construct+0x183>
   11fe3:	48 89 c7             	mov    %rax,%rdi
   11fe6:	e8 65 81 ff ff       	call   a150 <__assert_fail@plt>
   11feb:	48 83 7d b8 00       	cmpq   $0x0,-0x48(%rbp)
   11ff0:	75 28                	jne    1201a <_Z15alloc_particlesPP10particle_tiii+0x1cb>
   11ff2:	48 8d 05 ef 49 00 00 	lea    0x49ef(%rip),%rax        # 169e8 <_ZStL19piecewise_construct+0x288>
   11ff9:	48 89 c1             	mov    %rax,%rcx
   11ffc:	ba fc 00 00 00       	mov    $0xfc,%edx
   12001:	48 8d 05 98 47 00 00 	lea    0x4798(%rip),%rax        # 167a0 <_ZStL19piecewise_construct+0x40>
   12008:	48 89 c6             	mov    %rax,%rsi
   1200b:	48 8d 05 34 4a 00 00 	lea    0x4a34(%rip),%rax        # 16a46 <_ZStL19piecewise_construct+0x2e6>
   12012:	48 89 c7             	mov    %rax,%rdi
   12015:	e8 36 81 ff ff       	call   a150 <__assert_fail@plt>

    #pragma omp parallel for
   1201a:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   1201e:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   12022:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   12026:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   1202a:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   1202e:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   12032:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   12036:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   1203a:	8b 45 8c             	mov    -0x74(%rbp),%eax
   1203d:	89 45 f0             	mov    %eax,-0x10(%rbp)
   12040:	8b 45 90             	mov    -0x70(%rbp),%eax
   12043:	89 45 ec             	mov    %eax,-0x14(%rbp)
   12046:	8b 45 94             	mov    -0x6c(%rbp),%eax
   12049:	89 45 e8             	mov    %eax,-0x18(%rbp)
   1204c:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   12050:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   12054:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   12058:	b9 00 00 00 00       	mov    $0x0,%ecx
   1205d:	ba 00 00 00 00       	mov    $0x0,%edx
   12062:	48 89 c6             	mov    %rax,%rsi
   12065:	48 8d 05 ff 04 00 00 	lea    0x4ff(%rip),%rax        # 1256b <_Z15alloc_particlesPP10particle_tiii._omp_fn.0>
   1206c:	48 89 c7             	mov    %rax,%rdi
   1206f:	e8 3c 87 ff ff       	call   a7b0 <GOMP_parallel@plt>
        ppt->enabled = false;
        ppt->ndim = ndim;
        ppt->nfeat = nfeat;
        ppt->updateAcceleration = nullptr;
    }
}
   12074:	90                   	nop
   12075:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12079:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   12080:	00 00 
   12082:	74 05                	je     12089 <_Z15alloc_particlesPP10particle_tiii+0x23a>
   12084:	e8 57 84 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   12089:	c9                   	leave  
   1208a:	c3                   	ret    

000000000001208b <_Z15alloc_particlesP17chunk_particles_tiii>:

void alloc_particles(chunk_particles_t* chunkParticles, int nParticle, int ndim, int nfeat) {
   1208b:	f3 0f 1e fa          	endbr64 
   1208f:	55                   	push   %rbp
   12090:	48 89 e5             	mov    %rsp,%rbp
   12093:	48 83 ec 30          	sub    $0x30,%rsp
   12097:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   1209b:	89 75 e4             	mov    %esi,-0x1c(%rbp)
   1209e:	89 55 e0             	mov    %edx,-0x20(%rbp)
   120a1:	89 4d dc             	mov    %ecx,-0x24(%rbp)
   120a4:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   120ab:	00 00 
   120ad:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   120b1:	31 c0                	xor    %eax,%eax
    particle_t* particles = chunkParticles->particles;
   120b3:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   120b7:	48 8b 00             	mov    (%rax),%rax
   120ba:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
    alloc_particles(&particles, nParticle, ndim, nfeat);
   120be:	8b 4d dc             	mov    -0x24(%rbp),%ecx
   120c1:	8b 55 e0             	mov    -0x20(%rbp),%edx
   120c4:	8b 75 e4             	mov    -0x1c(%rbp),%esi
   120c7:	48 8d 45 f0          	lea    -0x10(%rbp),%rax
   120cb:	48 89 c7             	mov    %rax,%rdi
   120ce:	e8 7c fd ff ff       	call   11e4f <_Z15alloc_particlesPP10particle_tiii>
    chunkParticles->particles = particles;
   120d3:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   120d7:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   120db:	48 89 10             	mov    %rdx,(%rax)
    chunkParticles->nParticle = nParticle;
   120de:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   120e2:	8b 55 e4             	mov    -0x1c(%rbp),%edx
   120e5:	89 50 08             	mov    %edx,0x8(%rax)
}
   120e8:	90                   	nop
   120e9:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   120ed:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   120f4:	00 00 
   120f6:	74 05                	je     120fd <_Z15alloc_particlesP17chunk_particles_tiii+0x72>
   120f8:	e8 e3 83 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   120fd:	c9                   	leave  
   120fe:	c3                   	ret    

00000000000120ff <_Z41__static_initialization_and_destruction_0ii>:
   120ff:	f3 0f 1e fa          	endbr64 
   12103:	55                   	push   %rbp
   12104:	48 89 e5             	mov    %rsp,%rbp
   12107:	48 83 ec 10          	sub    $0x10,%rsp
   1210b:	89 7d fc             	mov    %edi,-0x4(%rbp)
   1210e:	89 75 f8             	mov    %esi,-0x8(%rbp)
   12111:	83 7d fc 01          	cmpl   $0x1,-0x4(%rbp)
   12115:	75 3b                	jne    12152 <_Z41__static_initialization_and_destruction_0ii+0x53>
   12117:	81 7d f8 ff ff 00 00 	cmpl   $0xffff,-0x8(%rbp)
   1211e:	75 32                	jne    12152 <_Z41__static_initialization_and_destruction_0ii+0x53>
   12120:	48 8d 05 2b c2 00 00 	lea    0xc22b(%rip),%rax        # 1e352 <_ZStL8__ioinit>
   12127:	48 89 c7             	mov    %rax,%rdi
   1212a:	e8 a1 86 ff ff       	call   a7d0 <_ZNSt8ios_base4InitC1Ev@plt>
   1212f:	48 8d 05 d2 ae 00 00 	lea    0xaed2(%rip),%rax        # 1d008 <__dso_handle>
   12136:	48 89 c2             	mov    %rax,%rdx
   12139:	48 8d 05 12 c2 00 00 	lea    0xc212(%rip),%rax        # 1e352 <_ZStL8__ioinit>
   12140:	48 89 c6             	mov    %rax,%rsi
   12143:	48 8b 05 ae ae 00 00 	mov    0xaeae(%rip),%rax        # 1cff8 <_ZNSt8ios_base4InitD1Ev@GLIBCXX_3.4>
   1214a:	48 89 c7             	mov    %rax,%rdi
   1214d:	e8 de 81 ff ff       	call   a330 <__cxa_atexit@plt>
   12152:	90                   	nop
   12153:	c9                   	leave  
   12154:	c3                   	ret    

0000000000012155 <_GLOBAL__sub_I__Z14updateVelocityP10particle_td>:
   12155:	f3 0f 1e fa          	endbr64 
   12159:	55                   	push   %rbp
   1215a:	48 89 e5             	mov    %rsp,%rbp
   1215d:	be ff ff 00 00       	mov    $0xffff,%esi
   12162:	bf 01 00 00 00       	mov    $0x1,%edi
   12167:	e8 93 ff ff ff       	call   120ff <_Z41__static_initialization_and_destruction_0ii>
   1216c:	5d                   	pop    %rbp
   1216d:	c3                   	ret    

000000000001216e <_Z18clearAccelerationsP10particle_ti._omp_fn.0>:
    #pragma omp parallel for
   1216e:	f3 0f 1e fa          	endbr64 
   12172:	55                   	push   %rbp
   12173:	48 89 e5             	mov    %rsp,%rbp
   12176:	41 54                	push   %r12
   12178:	53                   	push   %rbx
   12179:	48 83 ec 30          	sub    $0x30,%rsp
   1217d:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
   12181:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   12185:	8b 40 08             	mov    0x8(%rax),%eax
   12188:	89 45 e4             	mov    %eax,-0x1c(%rbp)
    for(int n=0; n<nParticle; n++) {
   1218b:	8b 5d e4             	mov    -0x1c(%rbp),%ebx
   1218e:	e8 ad 88 ff ff       	call   aa40 <omp_get_num_threads@plt>
   12193:	41 89 c4             	mov    %eax,%r12d
   12196:	e8 e5 86 ff ff       	call   a880 <omp_get_thread_num@plt>
   1219b:	89 c6                	mov    %eax,%esi
   1219d:	89 d8                	mov    %ebx,%eax
   1219f:	99                   	cltd   
   121a0:	41 f7 fc             	idiv   %r12d
   121a3:	89 c1                	mov    %eax,%ecx
   121a5:	89 d8                	mov    %ebx,%eax
   121a7:	99                   	cltd   
   121a8:	41 f7 fc             	idiv   %r12d
   121ab:	89 d0                	mov    %edx,%eax
   121ad:	39 c6                	cmp    %eax,%esi
   121af:	7c 6e                	jl     1221f <_Z18clearAccelerationsP10particle_ti._omp_fn.0+0xb1>
   121b1:	0f af f1             	imul   %ecx,%esi
   121b4:	89 f2                	mov    %esi,%edx
   121b6:	01 d0                	add    %edx,%eax
   121b8:	8d 14 08             	lea    (%rax,%rcx,1),%edx
   121bb:	39 d0                	cmp    %edx,%eax
   121bd:	7d 6a                	jge    12229 <_Z18clearAccelerationsP10particle_ti._omp_fn.0+0xbb>
   121bf:	89 45 dc             	mov    %eax,-0x24(%rbp)
        double* acceleration = pt->acceleration;
   121c2:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   121c6:	48 8b 00             	mov    (%rax),%rax
   121c9:	48 8b 40 18          	mov    0x18(%rax),%rax
   121cd:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
        for(int dim=0; dim<pt->ndim; dim++)
   121d1:	c7 45 e0 00 00 00 00 	movl   $0x0,-0x20(%rbp)
   121d8:	90                   	nop
   121d9:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   121dd:	48 8b 00             	mov    (%rax),%rax
   121e0:	8b 40 30             	mov    0x30(%rax),%eax
   121e3:	3b 45 e0             	cmp    -0x20(%rbp),%eax
   121e6:	7f 1d                	jg     12205 <_Z18clearAccelerationsP10particle_ti._omp_fn.0+0x97>
        pt++;
   121e8:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   121ec:	48 8b 00             	mov    (%rax),%rax
   121ef:	48 8d 48 40          	lea    0x40(%rax),%rcx
   121f3:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   121f7:	48 89 08             	mov    %rcx,(%rax)
   121fa:	83 45 dc 01          	addl   $0x1,-0x24(%rbp)
   121fe:	39 55 dc             	cmp    %edx,-0x24(%rbp)
   12201:	7c bf                	jl     121c2 <_Z18clearAccelerationsP10particle_ti._omp_fn.0+0x54>
    #pragma omp parallel for
   12203:	eb 24                	jmp    12229 <_Z18clearAccelerationsP10particle_ti._omp_fn.0+0xbb>
            *(acceleration++) = 0;
   12205:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12209:	48 8d 48 08          	lea    0x8(%rax),%rcx
   1220d:	48 89 4d e8          	mov    %rcx,-0x18(%rbp)
   12211:	66 0f ef c0          	pxor   %xmm0,%xmm0
   12215:	f2 0f 11 00          	movsd  %xmm0,(%rax)
        for(int dim=0; dim<pt->ndim; dim++)
   12219:	83 45 e0 01          	addl   $0x1,-0x20(%rbp)
   1221d:	eb ba                	jmp    121d9 <_Z18clearAccelerationsP10particle_ti._omp_fn.0+0x6b>
   1221f:	b8 00 00 00 00       	mov    $0x0,%eax
   12224:	83 c1 01             	add    $0x1,%ecx
   12227:	eb 88                	jmp    121b1 <_Z18clearAccelerationsP10particle_ti._omp_fn.0+0x43>
    #pragma omp parallel for
   12229:	90                   	nop
   1222a:	48 83 c4 30          	add    $0x30,%rsp
   1222e:	5b                   	pop    %rbx
   1222f:	41 5c                	pop    %r12
   12231:	5d                   	pop    %rbp
   12232:	c3                   	ret    

0000000000012233 <_Z16update_particlesP10particle_tid._omp_fn.0>:
    #pragma omp parallel for collapse(2)
   12233:	f3 0f 1e fa          	endbr64 
   12237:	55                   	push   %rbp
   12238:	48 89 e5             	mov    %rsp,%rbp
   1223b:	41 55                	push   %r13
   1223d:	41 54                	push   %r12
   1223f:	53                   	push   %rbx
   12240:	48 83 ec 38          	sub    $0x38,%rsp
   12244:	48 89 7d b8          	mov    %rdi,-0x48(%rbp)
   12248:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   1224c:	8b 40 08             	mov    0x8(%rax),%eax
   1224f:	89 45 d4             	mov    %eax,-0x2c(%rbp)
   12252:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   12256:	48 8b 00             	mov    (%rax),%rax
   12259:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
    for(int i=0; i<nParticle; i++) 
   1225d:	8b 55 d4             	mov    -0x2c(%rbp),%edx
   12260:	85 d2                	test   %edx,%edx
   12262:	0f 8e 23 02 00 00    	jle    1248b <_Z16update_particlesP10particle_tid._omp_fn.0+0x258>
   12268:	bb 00 00 00 00       	mov    $0x0,%ebx
   1226d:	41 bc 00 00 00 00    	mov    $0x0,%r12d
   12273:	89 d0                	mov    %edx,%eax
   12275:	44 29 e0             	sub    %r12d,%eax
   12278:	8d 48 ff             	lea    -0x1(%rax),%ecx
   1227b:	44 01 e1             	add    %r12d,%ecx
   1227e:	45 85 e4             	test   %r12d,%r12d
   12281:	0f 8f df 01 00 00    	jg     12466 <_Z16update_particlesP10particle_tid._omp_fn.0+0x233>
   12287:	e9 c7 01 00 00       	jmp    12453 <_Z16update_particlesP10particle_tid._omp_fn.0+0x220>
   1228c:	e8 af 87 ff ff       	call   aa40 <omp_get_num_threads@plt>
   12291:	41 89 c5             	mov    %eax,%r13d
   12294:	e8 e7 85 ff ff       	call   a880 <omp_get_thread_num@plt>
   12299:	89 c6                	mov    %eax,%esi
   1229b:	89 d8                	mov    %ebx,%eax
   1229d:	ba 00 00 00 00       	mov    $0x0,%edx
   122a2:	41 f7 f5             	div    %r13d
   122a5:	89 c1                	mov    %eax,%ecx
   122a7:	89 d8                	mov    %ebx,%eax
   122a9:	ba 00 00 00 00       	mov    $0x0,%edx
   122ae:	41 f7 f5             	div    %r13d
   122b1:	89 d0                	mov    %edx,%eax
   122b3:	39 c6                	cmp    %eax,%esi
   122b5:	0f 82 8b 01 00 00    	jb     12446 <_Z16update_particlesP10particle_tid._omp_fn.0+0x213>
   122bb:	0f af f1             	imul   %ecx,%esi
   122be:	89 f2                	mov    %esi,%edx
   122c0:	8d 1c 02             	lea    (%rdx,%rax,1),%ebx
   122c3:	44 8d 2c 0b          	lea    (%rbx,%rcx,1),%r13d
   122c7:	44 39 eb             	cmp    %r13d,%ebx
   122ca:	0f 83 bb 01 00 00    	jae    1248b <_Z16update_particlesP10particle_tid._omp_fn.0+0x258>
   122d0:	bf 00 00 00 00       	mov    $0x0,%edi
   122d5:	89 da                	mov    %ebx,%edx
   122d7:	49 63 c4             	movslq %r12d,%rax
   122da:	66 0f ef c9          	pxor   %xmm1,%xmm1
   122de:	f2 48 0f 2a c8       	cvtsi2sd %rax,%xmm1
   122e3:	48 85 d2             	test   %rdx,%rdx
   122e6:	78 0b                	js     122f3 <_Z16update_particlesP10particle_tid._omp_fn.0+0xc0>
   122e8:	66 0f ef c0          	pxor   %xmm0,%xmm0
   122ec:	f2 48 0f 2a c2       	cvtsi2sd %rdx,%xmm0
   122f1:	eb 1c                	jmp    1230f <_Z16update_particlesP10particle_tid._omp_fn.0+0xdc>
   122f3:	48 89 d0             	mov    %rdx,%rax
   122f6:	48 d1 e8             	shr    %rax
   122f9:	48 89 d1             	mov    %rdx,%rcx
   122fc:	83 e1 01             	and    $0x1,%ecx
   122ff:	48 09 c8             	or     %rcx,%rax
   12302:	66 0f ef c0          	pxor   %xmm0,%xmm0
   12306:	f2 48 0f 2a c0       	cvtsi2sd %rax,%xmm0
   1230b:	f2 0f 58 c0          	addsd  %xmm0,%xmm0
   1230f:	f2 0f 10 15 49 47 00 	movsd  0x4749(%rip),%xmm2        # 16a60 <_ZStL19piecewise_construct+0x300>
   12316:	00 
   12317:	f2 0f 5c ca          	subsd  %xmm2,%xmm1
   1231b:	66 0f 28 d0          	movapd %xmm0,%xmm2
   1231f:	f2 0f 58 d0          	addsd  %xmm0,%xmm2
   12323:	66 0f 28 c1          	movapd %xmm1,%xmm0
   12327:	f2 0f 59 c1          	mulsd  %xmm1,%xmm0
   1232b:	f2 0f 58 c2          	addsd  %xmm2,%xmm0
   1232f:	66 0f ef d2          	pxor   %xmm2,%xmm2
   12333:	66 0f 2f d0          	comisd %xmm0,%xmm2
   12337:	77 66                	ja     1239f <_Z16update_particlesP10particle_tid._omp_fn.0+0x16c>
   12339:	e9 8a 00 00 00       	jmp    123c8 <_Z16update_particlesP10particle_tid._omp_fn.0+0x195>
   1233e:	44 8b 65 cc          	mov    -0x34(%rbp),%r12d
            updateAcceleration(particles+i, particles+j);
   12342:	8b 45 d0             	mov    -0x30(%rbp),%eax
   12345:	48 98                	cltq   
   12347:	48 c1 e0 06          	shl    $0x6,%rax
   1234b:	48 89 c2             	mov    %rax,%rdx
   1234e:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12352:	48 01 c2             	add    %rax,%rdx
   12355:	8b 45 cc             	mov    -0x34(%rbp),%eax
   12358:	48 98                	cltq   
   1235a:	48 c1 e0 06          	shl    $0x6,%rax
   1235e:	48 89 c1             	mov    %rax,%rcx
   12361:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12365:	48 01 c8             	add    %rcx,%rax
   12368:	48 89 d6             	mov    %rdx,%rsi
   1236b:	48 89 c7             	mov    %rax,%rdi
   1236e:	e8 8c e7 ff ff       	call   10aff <_Z18updateAccelerationP10particle_tS0_>
   12373:	83 c3 01             	add    $0x1,%ebx
   12376:	44 39 eb             	cmp    %r13d,%ebx
   12379:	0f 83 0c 01 00 00    	jae    1248b <_Z16update_particlesP10particle_tid._omp_fn.0+0x258>
   1237f:	83 45 d0 01          	addl   $0x1,-0x30(%rbp)
   12383:	44 39 65 d0          	cmp    %r12d,-0x30(%rbp)
   12387:	7c b9                	jl     12342 <_Z16update_particlesP10particle_tid._omp_fn.0+0x10f>
   12389:	c7 45 d0 00 00 00 00 	movl   $0x0,-0x30(%rbp)
   12390:	83 45 cc 01          	addl   $0x1,-0x34(%rbp)
   12394:	44 8b 65 cc          	mov    -0x34(%rbp),%r12d
   12398:	45 85 e4             	test   %r12d,%r12d
   1239b:	7f a5                	jg     12342 <_Z16update_particlesP10particle_tid._omp_fn.0+0x10f>
   1239d:	eb ea                	jmp    12389 <_Z16update_particlesP10particle_tid._omp_fn.0+0x156>
   1239f:	b8 00 00 00 00       	mov    $0x0,%eax
   123a4:	85 c0                	test   %eax,%eax
   123a6:	7f 05                	jg     123ad <_Z16update_particlesP10particle_tid._omp_fn.0+0x17a>
   123a8:	83 c0 01             	add    $0x1,%eax
   123ab:	eb f7                	jmp    123a4 <_Z16update_particlesP10particle_tid._omp_fn.0+0x171>
   123ad:	89 c2                	mov    %eax,%edx
   123af:	01 fa                	add    %edi,%edx
   123b1:	39 da                	cmp    %ebx,%edx
   123b3:	77 04                	ja     123b9 <_Z16update_particlesP10particle_tid._omp_fn.0+0x186>
   123b5:	89 d7                	mov    %edx,%edi
   123b7:	eb ef                	jmp    123a8 <_Z16update_particlesP10particle_tid._omp_fn.0+0x175>
   123b9:	89 da                	mov    %ebx,%edx
   123bb:	29 fa                	sub    %edi,%edx
   123bd:	89 55 d0             	mov    %edx,-0x30(%rbp)
   123c0:	89 45 cc             	mov    %eax,-0x34(%rbp)
   123c3:	e9 76 ff ff ff       	jmp    1233e <_Z16update_particlesP10particle_tid._omp_fn.0+0x10b>
   123c8:	f2 0f 51 c0          	sqrtsd %xmm0,%xmm0
   123cc:	f2 0f 5c c1          	subsd  %xmm1,%xmm0
   123d0:	66 0f 2f 05 80 46 00 	comisd 0x4680(%rip),%xmm0        # 16a58 <_ZStL19piecewise_construct+0x2f8>
   123d7:	00 
   123d8:	73 07                	jae    123e1 <_Z16update_particlesP10particle_tid._omp_fn.0+0x1ae>
   123da:	f2 48 0f 2c c0       	cvttsd2si %xmm0,%rax
   123df:	eb 1e                	jmp    123ff <_Z16update_particlesP10particle_tid._omp_fn.0+0x1cc>
   123e1:	f2 0f 10 0d 6f 46 00 	movsd  0x466f(%rip),%xmm1        # 16a58 <_ZStL19piecewise_construct+0x2f8>
   123e8:	00 
   123e9:	f2 0f 5c c1          	subsd  %xmm1,%xmm0
   123ed:	f2 48 0f 2c c0       	cvttsd2si %xmm0,%rax
   123f2:	48 b9 00 00 00 00 00 	movabs $0x8000000000000000,%rcx
   123f9:	00 00 80 
   123fc:	48 31 c8             	xor    %rcx,%rax
   123ff:	48 8d 48 ff          	lea    -0x1(%rax),%rcx
   12403:	48 0f af c8          	imul   %rax,%rcx
   12407:	48 89 ce             	mov    %rcx,%rsi
   1240a:	48 d1 ee             	shr    %rsi
   1240d:	49 63 cc             	movslq %r12d,%rcx
   12410:	48 0f af c8          	imul   %rax,%rcx
   12414:	48 01 f1             	add    %rsi,%rcx
   12417:	49 63 f4             	movslq %r12d,%rsi
   1241a:	48 01 c6             	add    %rax,%rsi
   1241d:	48 39 ca             	cmp    %rcx,%rdx
   12420:	0f 82 79 ff ff ff    	jb     1239f <_Z16update_particlesP10particle_tid._omp_fn.0+0x16c>
   12426:	48 01 ce             	add    %rcx,%rsi
   12429:	48 39 f2             	cmp    %rsi,%rdx
   1242c:	0f 83 6d ff ff ff    	jae    1239f <_Z16update_particlesP10particle_tid._omp_fn.0+0x16c>
   12432:	44 01 e0             	add    %r12d,%eax
   12435:	89 45 cc             	mov    %eax,-0x34(%rbp)
   12438:	48 89 d0             	mov    %rdx,%rax
   1243b:	48 29 c8             	sub    %rcx,%rax
   1243e:	89 45 d0             	mov    %eax,-0x30(%rbp)
   12441:	e9 f8 fe ff ff       	jmp    1233e <_Z16update_particlesP10particle_tid._omp_fn.0+0x10b>
   12446:	b8 00 00 00 00       	mov    $0x0,%eax
   1244b:	83 c1 01             	add    $0x1,%ecx
   1244e:	e9 68 fe ff ff       	jmp    122bb <_Z16update_particlesP10particle_tid._omp_fn.0+0x88>
   12453:	85 c9                	test   %ecx,%ecx
   12455:	0f 8e 31 fe ff ff    	jle    1228c <_Z16update_particlesP10particle_tid._omp_fn.0+0x59>
   1245b:	41 bc 01 00 00 00    	mov    $0x1,%r12d
   12461:	e9 0d fe ff ff       	jmp    12273 <_Z16update_particlesP10particle_tid._omp_fn.0+0x40>
   12466:	85 c9                	test   %ecx,%ecx
   12468:	7f 0a                	jg     12474 <_Z16update_particlesP10particle_tid._omp_fn.0+0x241>
   1246a:	ba 00 00 00 00       	mov    $0x0,%edx
   1246f:	e9 ff fd ff ff       	jmp    12273 <_Z16update_particlesP10particle_tid._omp_fn.0+0x40>
   12474:	89 c1                	mov    %eax,%ecx
   12476:	41 0f af cc          	imul   %r12d,%ecx
   1247a:	8d 50 ff             	lea    -0x1(%rax),%edx
   1247d:	0f af c2             	imul   %edx,%eax
   12480:	d1 f8                	sar    %eax
   12482:	01 c8                	add    %ecx,%eax
   12484:	89 c3                	mov    %eax,%ebx
   12486:	e9 01 fe ff ff       	jmp    1228c <_Z16update_particlesP10particle_tid._omp_fn.0+0x59>
    #pragma omp parallel for collapse(2)
   1248b:	90                   	nop
   1248c:	48 83 c4 38          	add    $0x38,%rsp
   12490:	5b                   	pop    %rbx
   12491:	41 5c                	pop    %r12
   12493:	41 5d                	pop    %r13
   12495:	5d                   	pop    %rbp
   12496:	c3                   	ret    

0000000000012497 <_Z16update_particlesP10particle_tid._omp_fn.1>:
    #pragma omp parallel for
   12497:	f3 0f 1e fa          	endbr64 
   1249b:	55                   	push   %rbp
   1249c:	48 89 e5             	mov    %rsp,%rbp
   1249f:	41 54                	push   %r12
   124a1:	53                   	push   %rbx
   124a2:	48 83 ec 30          	sub    $0x30,%rsp
   124a6:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
   124aa:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   124ae:	f2 0f 10 40 08       	movsd  0x8(%rax),%xmm0
   124b3:	f2 0f 11 45 e0       	movsd  %xmm0,-0x20(%rbp)
   124b8:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   124bc:	8b 40 10             	mov    0x10(%rax),%eax
   124bf:	89 45 dc             	mov    %eax,-0x24(%rbp)
   124c2:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   124c6:	48 8b 00             	mov    (%rax),%rax
   124c9:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    for(int i=0; i<nParticle; i++) {
   124cd:	8b 5d dc             	mov    -0x24(%rbp),%ebx
   124d0:	e8 6b 85 ff ff       	call   aa40 <omp_get_num_threads@plt>
   124d5:	41 89 c4             	mov    %eax,%r12d
   124d8:	e8 a3 83 ff ff       	call   a880 <omp_get_thread_num@plt>
   124dd:	89 c6                	mov    %eax,%esi
   124df:	89 d8                	mov    %ebx,%eax
   124e1:	99                   	cltd   
   124e2:	41 f7 fc             	idiv   %r12d
   124e5:	89 c1                	mov    %eax,%ecx
   124e7:	89 d8                	mov    %ebx,%eax
   124e9:	99                   	cltd   
   124ea:	41 f7 fc             	idiv   %r12d
   124ed:	89 d0                	mov    %edx,%eax
   124ef:	39 c6                	cmp    %eax,%esi
   124f1:	7c 64                	jl     12557 <_Z16update_particlesP10particle_tid._omp_fn.1+0xc0>
   124f3:	0f af f1             	imul   %ecx,%esi
   124f6:	89 f2                	mov    %esi,%edx
   124f8:	01 d0                	add    %edx,%eax
   124fa:	8d 1c 08             	lea    (%rax,%rcx,1),%ebx
   124fd:	39 d8                	cmp    %ebx,%eax
   124ff:	7d 60                	jge    12561 <_Z16update_particlesP10particle_tid._omp_fn.1+0xca>
   12501:	89 45 d8             	mov    %eax,-0x28(%rbp)
        updateVelocity(particles+i, timeStep);
   12504:	8b 45 d8             	mov    -0x28(%rbp),%eax
   12507:	48 98                	cltq   
   12509:	48 c1 e0 06          	shl    $0x6,%rax
   1250d:	48 89 c2             	mov    %rax,%rdx
   12510:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12514:	48 01 c2             	add    %rax,%rdx
   12517:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   1251b:	66 48 0f 6e c0       	movq   %rax,%xmm0
   12520:	48 89 d7             	mov    %rdx,%rdi
   12523:	e8 c9 e4 ff ff       	call   109f1 <_Z14updateVelocityP10particle_td>
        updatePosition(particles+i, timeStep);
   12528:	8b 45 d8             	mov    -0x28(%rbp),%eax
   1252b:	48 98                	cltq   
   1252d:	48 c1 e0 06          	shl    $0x6,%rax
   12531:	48 89 c2             	mov    %rax,%rdx
   12534:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12538:	48 01 c2             	add    %rax,%rdx
   1253b:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   1253f:	66 48 0f 6e c0       	movq   %rax,%xmm0
   12544:	48 89 d7             	mov    %rdx,%rdi
   12547:	e8 2c e5 ff ff       	call   10a78 <_Z14updatePositionP10particle_td>
   1254c:	83 45 d8 01          	addl   $0x1,-0x28(%rbp)
   12550:	39 5d d8             	cmp    %ebx,-0x28(%rbp)
   12553:	7c af                	jl     12504 <_Z16update_particlesP10particle_tid._omp_fn.1+0x6d>
    #pragma omp parallel for
   12555:	eb 0a                	jmp    12561 <_Z16update_particlesP10particle_tid._omp_fn.1+0xca>
   12557:	b8 00 00 00 00       	mov    $0x0,%eax
   1255c:	83 c1 01             	add    $0x1,%ecx
   1255f:	eb 92                	jmp    124f3 <_Z16update_particlesP10particle_tid._omp_fn.1+0x5c>
   12561:	90                   	nop
   12562:	48 83 c4 30          	add    $0x30,%rsp
   12566:	5b                   	pop    %rbx
   12567:	41 5c                	pop    %r12
   12569:	5d                   	pop    %rbp
   1256a:	c3                   	ret    

000000000001256b <_Z15alloc_particlesPP10particle_tiii._omp_fn.0>:
    #pragma omp parallel for
   1256b:	f3 0f 1e fa          	endbr64 
   1256f:	55                   	push   %rbp
   12570:	48 89 e5             	mov    %rsp,%rbp
   12573:	41 54                	push   %r12
   12575:	53                   	push   %rbx
   12576:	48 83 ec 50          	sub    $0x50,%rsp
   1257a:	48 89 7d a8          	mov    %rdi,-0x58(%rbp)
   1257e:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12582:	48 8b 40 20          	mov    0x20(%rax),%rax
   12586:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   1258a:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   1258e:	48 8b 40 18          	mov    0x18(%rax),%rax
   12592:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   12596:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   1259a:	48 8b 40 10          	mov    0x10(%rax),%rax
   1259e:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   125a2:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   125a6:	48 8b 40 08          	mov    0x8(%rax),%rax
   125aa:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   125ae:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   125b2:	8b 40 30             	mov    0x30(%rax),%eax
   125b5:	89 45 b4             	mov    %eax,-0x4c(%rbp)
   125b8:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   125bc:	8b 40 2c             	mov    0x2c(%rax),%eax
   125bf:	89 45 b8             	mov    %eax,-0x48(%rbp)
   125c2:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   125c6:	8b 40 28             	mov    0x28(%rax),%eax
   125c9:	89 45 bc             	mov    %eax,-0x44(%rbp)
   125cc:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   125d0:	48 8b 00             	mov    (%rax),%rax
   125d3:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
    for(int n=0; n<nParticle; n++) {
   125d7:	8b 5d bc             	mov    -0x44(%rbp),%ebx
   125da:	e8 61 84 ff ff       	call   aa40 <omp_get_num_threads@plt>
   125df:	41 89 c4             	mov    %eax,%r12d
   125e2:	e8 99 82 ff ff       	call   a880 <omp_get_thread_num@plt>
   125e7:	89 c6                	mov    %eax,%esi
   125e9:	89 d8                	mov    %ebx,%eax
   125eb:	99                   	cltd   
   125ec:	41 f7 fc             	idiv   %r12d
   125ef:	89 c1                	mov    %eax,%ecx
   125f1:	89 d8                	mov    %ebx,%eax
   125f3:	99                   	cltd   
   125f4:	41 f7 fc             	idiv   %r12d
   125f7:	89 d0                	mov    %edx,%eax
   125f9:	39 c6                	cmp    %eax,%esi
   125fb:	0f 8c ee 00 00 00    	jl     126ef <_Z15alloc_particlesPP10particle_tiii._omp_fn.0+0x184>
   12601:	0f af f1             	imul   %ecx,%esi
   12604:	89 f2                	mov    %esi,%edx
   12606:	01 d0                	add    %edx,%eax
   12608:	8d 14 08             	lea    (%rax,%rcx,1),%edx
   1260b:	39 d0                	cmp    %edx,%eax
   1260d:	0f 8d e9 00 00 00    	jge    126fc <_Z15alloc_particlesPP10particle_tiii._omp_fn.0+0x191>
   12613:	89 45 b0             	mov    %eax,-0x50(%rbp)
        particle_t* ppt = (*particles) + n;
   12616:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   1261a:	48 8b 08             	mov    (%rax),%rcx
   1261d:	8b 45 b0             	mov    -0x50(%rbp),%eax
   12620:	48 98                	cltq   
   12622:	48 c1 e0 06          	shl    $0x6,%rax
   12626:	48 01 c8             	add    %rcx,%rax
   12629:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
        ppt->position = pos + n*ndim;
   1262d:	8b 45 b0             	mov    -0x50(%rbp),%eax
   12630:	0f af 45 b8          	imul   -0x48(%rbp),%eax
   12634:	48 98                	cltq   
   12636:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   1263d:	00 
   1263e:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12642:	48 01 c1             	add    %rax,%rcx
   12645:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12649:	48 89 48 08          	mov    %rcx,0x8(%rax)
        ppt->velocity = vel + n*ndim;
   1264d:	8b 45 b0             	mov    -0x50(%rbp),%eax
   12650:	0f af 45 b8          	imul   -0x48(%rbp),%eax
   12654:	48 98                	cltq   
   12656:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   1265d:	00 
   1265e:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   12662:	48 01 c1             	add    %rax,%rcx
   12665:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12669:	48 89 48 10          	mov    %rcx,0x10(%rax)
        ppt->acceleration = acc + n*ndim;
   1266d:	8b 45 b0             	mov    -0x50(%rbp),%eax
   12670:	0f af 45 b8          	imul   -0x48(%rbp),%eax
   12674:	48 98                	cltq   
   12676:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   1267d:	00 
   1267e:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   12682:	48 01 c1             	add    %rax,%rcx
   12685:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12689:	48 89 48 18          	mov    %rcx,0x18(%rax)
        ppt->features = feat + n*nfeat;
   1268d:	8b 45 b0             	mov    -0x50(%rbp),%eax
   12690:	0f af 45 b4          	imul   -0x4c(%rbp),%eax
   12694:	48 98                	cltq   
   12696:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   1269d:	00 
   1269e:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   126a2:	48 01 c1             	add    %rax,%rcx
   126a5:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   126a9:	48 89 48 20          	mov    %rcx,0x20(%rax)
        ppt->id = reinterpret_cast<uint64_t>(ppt);
   126ad:	48 8b 4d e8          	mov    -0x18(%rbp),%rcx
   126b1:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   126b5:	48 89 08             	mov    %rcx,(%rax)
        ppt->enabled = false;
   126b8:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   126bc:	c6 40 38 00          	movb   $0x0,0x38(%rax)
        ppt->ndim = ndim;
   126c0:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   126c4:	8b 4d b8             	mov    -0x48(%rbp),%ecx
   126c7:	89 48 30             	mov    %ecx,0x30(%rax)
        ppt->nfeat = nfeat;
   126ca:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   126ce:	8b 4d b4             	mov    -0x4c(%rbp),%ecx
   126d1:	89 48 34             	mov    %ecx,0x34(%rax)
        ppt->updateAcceleration = nullptr;
   126d4:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   126d8:	48 c7 40 28 00 00 00 	movq   $0x0,0x28(%rax)
   126df:	00 
   126e0:	83 45 b0 01          	addl   $0x1,-0x50(%rbp)
   126e4:	39 55 b0             	cmp    %edx,-0x50(%rbp)
   126e7:	0f 8c 29 ff ff ff    	jl     12616 <_Z15alloc_particlesPP10particle_tiii._omp_fn.0+0xab>
    #pragma omp parallel for
   126ed:	eb 0d                	jmp    126fc <_Z15alloc_particlesPP10particle_tiii._omp_fn.0+0x191>
   126ef:	b8 00 00 00 00       	mov    $0x0,%eax
   126f4:	83 c1 01             	add    $0x1,%ecx
   126f7:	e9 05 ff ff ff       	jmp    12601 <_Z15alloc_particlesPP10particle_tiii._omp_fn.0+0x96>
   126fc:	90                   	nop
   126fd:	48 83 c4 50          	add    $0x50,%rsp
   12701:	5b                   	pop    %rbx
   12702:	41 5c                	pop    %r12
   12704:	5d                   	pop    %rbp
   12705:	c3                   	ret    

0000000000012706 <_ZStanSt13_Ios_FmtflagsS_>:
      _S_ios_fmtflags_min = ~__INT_MAX__
    };

  inline _GLIBCXX_CONSTEXPR _Ios_Fmtflags
  operator&(_Ios_Fmtflags __a, _Ios_Fmtflags __b)
  { return _Ios_Fmtflags(static_cast<int>(__a) & static_cast<int>(__b)); }
   12706:	f3 0f 1e fa          	endbr64 
   1270a:	55                   	push   %rbp
   1270b:	48 89 e5             	mov    %rsp,%rbp
   1270e:	89 7d fc             	mov    %edi,-0x4(%rbp)
   12711:	89 75 f8             	mov    %esi,-0x8(%rbp)
   12714:	8b 45 fc             	mov    -0x4(%rbp),%eax
   12717:	23 45 f8             	and    -0x8(%rbp),%eax
   1271a:	5d                   	pop    %rbp
   1271b:	c3                   	ret    

000000000001271c <_ZStorSt13_Ios_FmtflagsS_>:

  inline _GLIBCXX_CONSTEXPR _Ios_Fmtflags
  operator|(_Ios_Fmtflags __a, _Ios_Fmtflags __b)
  { return _Ios_Fmtflags(static_cast<int>(__a) | static_cast<int>(__b)); }
   1271c:	f3 0f 1e fa          	endbr64 
   12720:	55                   	push   %rbp
   12721:	48 89 e5             	mov    %rsp,%rbp
   12724:	89 7d fc             	mov    %edi,-0x4(%rbp)
   12727:	89 75 f8             	mov    %esi,-0x8(%rbp)
   1272a:	8b 45 fc             	mov    -0x4(%rbp),%eax
   1272d:	0b 45 f8             	or     -0x8(%rbp),%eax
   12730:	5d                   	pop    %rbp
   12731:	c3                   	ret    

0000000000012732 <_ZStcoSt13_Ios_Fmtflags>:
  operator^(_Ios_Fmtflags __a, _Ios_Fmtflags __b)
  { return _Ios_Fmtflags(static_cast<int>(__a) ^ static_cast<int>(__b)); }

  inline _GLIBCXX_CONSTEXPR _Ios_Fmtflags
  operator~(_Ios_Fmtflags __a)
  { return _Ios_Fmtflags(~static_cast<int>(__a)); }
   12732:	f3 0f 1e fa          	endbr64 
   12736:	55                   	push   %rbp
   12737:	48 89 e5             	mov    %rsp,%rbp
   1273a:	89 7d fc             	mov    %edi,-0x4(%rbp)
   1273d:	8b 45 fc             	mov    -0x4(%rbp),%eax
   12740:	f7 d0                	not    %eax
   12742:	5d                   	pop    %rbp
   12743:	c3                   	ret    

0000000000012744 <_ZStoRRSt13_Ios_FmtflagsS_>:

  inline const _Ios_Fmtflags&
  operator|=(_Ios_Fmtflags& __a, _Ios_Fmtflags __b)
  { return __a = __a | __b; }
   12744:	f3 0f 1e fa          	endbr64 
   12748:	55                   	push   %rbp
   12749:	48 89 e5             	mov    %rsp,%rbp
   1274c:	48 83 ec 10          	sub    $0x10,%rsp
   12750:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   12754:	89 75 f4             	mov    %esi,-0xc(%rbp)
   12757:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1275b:	8b 00                	mov    (%rax),%eax
   1275d:	8b 55 f4             	mov    -0xc(%rbp),%edx
   12760:	89 d6                	mov    %edx,%esi
   12762:	89 c7                	mov    %eax,%edi
   12764:	e8 b3 ff ff ff       	call   1271c <_ZStorSt13_Ios_FmtflagsS_>
   12769:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   1276d:	89 02                	mov    %eax,(%rdx)
   1276f:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12773:	c9                   	leave  
   12774:	c3                   	ret    

0000000000012775 <_ZStaNRSt13_Ios_FmtflagsS_>:

  inline const _Ios_Fmtflags&
  operator&=(_Ios_Fmtflags& __a, _Ios_Fmtflags __b)
  { return __a = __a & __b; }
   12775:	f3 0f 1e fa          	endbr64 
   12779:	55                   	push   %rbp
   1277a:	48 89 e5             	mov    %rsp,%rbp
   1277d:	48 83 ec 10          	sub    $0x10,%rsp
   12781:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   12785:	89 75 f4             	mov    %esi,-0xc(%rbp)
   12788:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1278c:	8b 00                	mov    (%rax),%eax
   1278e:	8b 55 f4             	mov    -0xc(%rbp),%edx
   12791:	89 d6                	mov    %edx,%esi
   12793:	89 c7                	mov    %eax,%edi
   12795:	e8 6c ff ff ff       	call   12706 <_ZStanSt13_Ios_FmtflagsS_>
   1279a:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   1279e:	89 02                	mov    %eax,(%rdx)
   127a0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   127a4:	c9                   	leave  
   127a5:	c3                   	ret    

00000000000127a6 <_ZStorSt13_Ios_OpenmodeS_>:
  operator&(_Ios_Openmode __a, _Ios_Openmode __b)
  { return _Ios_Openmode(static_cast<int>(__a) & static_cast<int>(__b)); }

  inline _GLIBCXX_CONSTEXPR _Ios_Openmode
  operator|(_Ios_Openmode __a, _Ios_Openmode __b)
  { return _Ios_Openmode(static_cast<int>(__a) | static_cast<int>(__b)); }
   127a6:	f3 0f 1e fa          	endbr64 
   127aa:	55                   	push   %rbp
   127ab:	48 89 e5             	mov    %rsp,%rbp
   127ae:	89 7d fc             	mov    %edi,-0x4(%rbp)
   127b1:	89 75 f8             	mov    %esi,-0x8(%rbp)
   127b4:	8b 45 fc             	mov    -0x4(%rbp),%eax
   127b7:	0b 45 f8             	or     -0x8(%rbp),%eax
   127ba:	5d                   	pop    %rbp
   127bb:	c3                   	ret    

00000000000127bc <_ZNSt8ios_base4setfESt13_Ios_FmtflagsS0_>:
     *
     *  This function clears @a mask in the format flags, then sets
     *  @a fmtfl @c & @a mask.  An example mask is @c ios_base::adjustfield.
    */
    fmtflags
    setf(fmtflags __fmtfl, fmtflags __mask)
   127bc:	f3 0f 1e fa          	endbr64 
   127c0:	55                   	push   %rbp
   127c1:	48 89 e5             	mov    %rsp,%rbp
   127c4:	48 83 ec 20          	sub    $0x20,%rsp
   127c8:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   127cc:	89 75 e4             	mov    %esi,-0x1c(%rbp)
   127cf:	89 55 e0             	mov    %edx,-0x20(%rbp)
    {
      fmtflags __old = _M_flags;
   127d2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   127d6:	8b 40 18             	mov    0x18(%rax),%eax
   127d9:	89 45 fc             	mov    %eax,-0x4(%rbp)
      _M_flags &= ~__mask;
   127dc:	8b 45 e0             	mov    -0x20(%rbp),%eax
   127df:	89 c7                	mov    %eax,%edi
   127e1:	e8 4c ff ff ff       	call   12732 <_ZStcoSt13_Ios_Fmtflags>
   127e6:	89 c2                	mov    %eax,%edx
   127e8:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   127ec:	48 83 c0 18          	add    $0x18,%rax
   127f0:	89 d6                	mov    %edx,%esi
   127f2:	48 89 c7             	mov    %rax,%rdi
   127f5:	e8 7b ff ff ff       	call   12775 <_ZStaNRSt13_Ios_FmtflagsS_>
      _M_flags |= (__fmtfl & __mask);
   127fa:	8b 55 e0             	mov    -0x20(%rbp),%edx
   127fd:	8b 45 e4             	mov    -0x1c(%rbp),%eax
   12800:	89 d6                	mov    %edx,%esi
   12802:	89 c7                	mov    %eax,%edi
   12804:	e8 fd fe ff ff       	call   12706 <_ZStanSt13_Ios_FmtflagsS_>
   12809:	89 c2                	mov    %eax,%edx
   1280b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1280f:	48 83 c0 18          	add    $0x18,%rax
   12813:	89 d6                	mov    %edx,%esi
   12815:	48 89 c7             	mov    %rax,%rdi
   12818:	e8 27 ff ff ff       	call   12744 <_ZStoRRSt13_Ios_FmtflagsS_>
      return __old;
   1281d:	8b 45 fc             	mov    -0x4(%rbp),%eax
    }
   12820:	c9                   	leave  
   12821:	c3                   	ret    

0000000000012822 <_ZSt10scientificRSt8ios_base>:
  }

  /// Calls base.setf(ios_base::scientific, ios_base::floatfield).
  inline ios_base&
  scientific(ios_base& __base)
  {
   12822:	f3 0f 1e fa          	endbr64 
   12826:	55                   	push   %rbp
   12827:	48 89 e5             	mov    %rsp,%rbp
   1282a:	48 83 ec 10          	sub    $0x10,%rsp
   1282e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    __base.setf(ios_base::scientific, ios_base::floatfield);
   12832:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12836:	ba 04 01 00 00       	mov    $0x104,%edx
   1283b:	be 00 01 00 00       	mov    $0x100,%esi
   12840:	48 89 c7             	mov    %rax,%rdi
   12843:	e8 74 ff ff ff       	call   127bc <_ZNSt8ios_base4setfESt13_Ios_FmtflagsS0_>
    return __base;
   12848:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
  }
   1284c:	c9                   	leave  
   1284d:	c3                   	ret    

000000000001284e <_Z18setupCommunicatorsP10topology_tii>:

__inline__ int sizeAfterPadding(int nParticles, int minimalGrid) {
    return (nParticles+minimalGrid-1)/minimalGrid * minimalGrid;
}

void setupCommunicators(topology_t* topology, int gridX, int gridY) {
   1284e:	f3 0f 1e fa          	endbr64 
   12852:	55                   	push   %rbp
   12853:	48 89 e5             	mov    %rsp,%rbp
   12856:	48 83 ec 20          	sub    $0x20,%rsp
   1285a:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   1285e:	89 75 e4             	mov    %esi,-0x1c(%rbp)
   12861:	89 55 e0             	mov    %edx,-0x20(%rbp)
   12864:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   1286b:	00 00 
   1286d:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   12871:	31 c0                	xor    %eax,%eax
    topology->gridX = gridX;
   12873:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12877:	8b 55 e4             	mov    -0x1c(%rbp),%edx
   1287a:	89 50 08             	mov    %edx,0x8(%rax)
    topology->gridY = gridY;
   1287d:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12881:	8b 55 e0             	mov    -0x20(%rbp),%edx
   12884:	89 50 0c             	mov    %edx,0xc(%rax)
    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   12887:	48 8d 45 f0          	lea    -0x10(%rbp),%rax
   1288b:	48 89 c6             	mov    %rax,%rsi
   1288e:	48 8d 05 8b b1 00 00 	lea    0xb18b(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   12895:	48 89 c7             	mov    %rax,%rdi
   12898:	e8 53 7c ff ff       	call   a4f0 <MPI_Comm_size@plt>
    assert(world_size == gridX*gridY);
   1289d:	8b 45 e4             	mov    -0x1c(%rbp),%eax
   128a0:	0f af 45 e0          	imul   -0x20(%rbp),%eax
   128a4:	89 c2                	mov    %eax,%edx
   128a6:	8b 45 f0             	mov    -0x10(%rbp),%eax
   128a9:	39 c2                	cmp    %eax,%edx
   128ab:	74 28                	je     128d5 <_Z18setupCommunicatorsP10topology_tii+0x87>
   128ad:	48 8d 05 f4 43 00 00 	lea    0x43f4(%rip),%rax        # 16ca8 <_ZN3MPIL11LOCK_SHAREDE+0x8>
   128b4:	48 89 c1             	mov    %rax,%rcx
   128b7:	ba 0f 00 00 00       	mov    $0xf,%edx
   128bc:	48 8d 05 15 44 00 00 	lea    0x4415(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   128c3:	48 89 c6             	mov    %rax,%rsi
   128c6:	48 8d 05 63 44 00 00 	lea    0x4463(%rip),%rax        # 16d30 <_ZN3MPIL11LOCK_SHAREDE+0x90>
   128cd:	48 89 c7             	mov    %rax,%rdi
   128d0:	e8 7b 78 ff ff       	call   a150 <__assert_fail@plt>
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   128d5:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
   128d9:	48 89 c6             	mov    %rax,%rsi
   128dc:	48 8d 05 3d b1 00 00 	lea    0xb13d(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   128e3:	48 89 c7             	mov    %rax,%rdi
   128e6:	e8 15 7e ff ff       	call   a700 <MPI_Comm_rank@plt>
    topology->rankX = rank % gridX;
   128eb:	8b 45 f4             	mov    -0xc(%rbp),%eax
   128ee:	99                   	cltd   
   128ef:	f7 7d e4             	idivl  -0x1c(%rbp)
   128f2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   128f6:	89 10                	mov    %edx,(%rax)
    topology->rankY = rank / gridX;
   128f8:	8b 45 f4             	mov    -0xc(%rbp),%eax
   128fb:	99                   	cltd   
   128fc:	f7 7d e4             	idivl  -0x1c(%rbp)
   128ff:	89 c2                	mov    %eax,%edx
   12901:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12905:	89 50 04             	mov    %edx,0x4(%rax)

    if(topology->rankY==0) 
   12908:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1290c:	8b 40 04             	mov    0x4(%rax),%eax
   1290f:	85 c0                	test   %eax,%eax
   12911:	75 28                	jne    1293b <_Z18setupCommunicatorsP10topology_tii+0xed>
        MPI_Comm_split(MPI_COMM_WORLD, topology->rankY, topology->rankX, &(topology->scatterXComm));
   12913:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12917:	48 8d 48 20          	lea    0x20(%rax),%rcx
   1291b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1291f:	8b 10                	mov    (%rax),%edx
   12921:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12925:	8b 40 04             	mov    0x4(%rax),%eax
   12928:	89 c6                	mov    %eax,%esi
   1292a:	48 8d 05 ef b0 00 00 	lea    0xb0ef(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   12931:	48 89 c7             	mov    %rax,%rdi
   12934:	e8 57 7f ff ff       	call   a890 <MPI_Comm_split@plt>
   12939:	eb 27                	jmp    12962 <_Z18setupCommunicatorsP10topology_tii+0x114>
    else
        MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, topology->rankX, &(topology->scatterXComm));
   1293b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1293f:	48 8d 50 20          	lea    0x20(%rax),%rdx
   12943:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12947:	8b 00                	mov    (%rax),%eax
   12949:	48 89 d1             	mov    %rdx,%rcx
   1294c:	89 c2                	mov    %eax,%edx
   1294e:	be 02 80 ff ff       	mov    $0xffff8002,%esi
   12953:	48 8d 05 c6 b0 00 00 	lea    0xb0c6(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   1295a:	48 89 c7             	mov    %rax,%rdi
   1295d:	e8 2e 7f ff ff       	call   a890 <MPI_Comm_split@plt>
    MPI_Barrier(MPI_COMM_WORLD);
   12962:	48 8d 05 b7 b0 00 00 	lea    0xb0b7(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   12969:	48 89 c7             	mov    %rax,%rdi
   1296c:	e8 6f 78 ff ff       	call   a1e0 <MPI_Barrier@plt>

    if(topology->rankX==0)
   12971:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12975:	8b 00                	mov    (%rax),%eax
   12977:	85 c0                	test   %eax,%eax
   12979:	75 28                	jne    129a3 <_Z18setupCommunicatorsP10topology_tii+0x155>
        MPI_Comm_split(MPI_COMM_WORLD, topology->rankX, topology->rankY, &(topology->scatterYComm));
   1297b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1297f:	48 8d 48 28          	lea    0x28(%rax),%rcx
   12983:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12987:	8b 50 04             	mov    0x4(%rax),%edx
   1298a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1298e:	8b 00                	mov    (%rax),%eax
   12990:	89 c6                	mov    %eax,%esi
   12992:	48 8d 05 87 b0 00 00 	lea    0xb087(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   12999:	48 89 c7             	mov    %rax,%rdi
   1299c:	e8 ef 7e ff ff       	call   a890 <MPI_Comm_split@plt>
   129a1:	eb 28                	jmp    129cb <_Z18setupCommunicatorsP10topology_tii+0x17d>
    else
        MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, topology->rankY, &(topology->scatterYComm));
   129a3:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   129a7:	48 8d 50 28          	lea    0x28(%rax),%rdx
   129ab:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   129af:	8b 40 04             	mov    0x4(%rax),%eax
   129b2:	48 89 d1             	mov    %rdx,%rcx
   129b5:	89 c2                	mov    %eax,%edx
   129b7:	be 02 80 ff ff       	mov    $0xffff8002,%esi
   129bc:	48 8d 05 5d b0 00 00 	lea    0xb05d(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   129c3:	48 89 c7             	mov    %rax,%rdi
   129c6:	e8 c5 7e ff ff       	call   a890 <MPI_Comm_split@plt>
    MPI_Barrier(MPI_COMM_WORLD);
   129cb:	48 8d 05 4e b0 00 00 	lea    0xb04e(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   129d2:	48 89 c7             	mov    %rax,%rdi
   129d5:	e8 06 78 ff ff       	call   a1e0 <MPI_Barrier@plt>

    MPI_Comm_split(MPI_COMM_WORLD, topology->rankX, topology->rankY, &(topology->reduceXComm));
   129da:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   129de:	48 8d 48 10          	lea    0x10(%rax),%rcx
   129e2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   129e6:	8b 50 04             	mov    0x4(%rax),%edx
   129e9:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   129ed:	8b 00                	mov    (%rax),%eax
   129ef:	89 c6                	mov    %eax,%esi
   129f1:	48 8d 05 28 b0 00 00 	lea    0xb028(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   129f8:	48 89 c7             	mov    %rax,%rdi
   129fb:	e8 90 7e ff ff       	call   a890 <MPI_Comm_split@plt>
    MPI_Barrier(MPI_COMM_WORLD);
   12a00:	48 8d 05 19 b0 00 00 	lea    0xb019(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   12a07:	48 89 c7             	mov    %rax,%rdi
   12a0a:	e8 d1 77 ff ff       	call   a1e0 <MPI_Barrier@plt>

    MPI_Comm_split(MPI_COMM_WORLD, topology->rankY, topology->rankX, &(topology->reduceYComm));
   12a0f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12a13:	48 8d 48 18          	lea    0x18(%rax),%rcx
   12a17:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12a1b:	8b 10                	mov    (%rax),%edx
   12a1d:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12a21:	8b 40 04             	mov    0x4(%rax),%eax
   12a24:	89 c6                	mov    %eax,%esi
   12a26:	48 8d 05 f3 af 00 00 	lea    0xaff3(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   12a2d:	48 89 c7             	mov    %rax,%rdi
   12a30:	e8 5b 7e ff ff       	call   a890 <MPI_Comm_split@plt>
    MPI_Barrier(MPI_COMM_WORLD);
   12a35:	48 8d 05 e4 af 00 00 	lea    0xafe4(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   12a3c:	48 89 c7             	mov    %rax,%rdi
   12a3f:	e8 9c 77 ff ff       	call   a1e0 <MPI_Barrier@plt>
}
   12a44:	90                   	nop
   12a45:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12a49:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   12a50:	00 00 
   12a52:	74 05                	je     12a59 <_Z18setupCommunicatorsP10topology_tii+0x20b>
   12a54:	e8 87 7a ff ff       	call   a4e0 <__stack_chk_fail@plt>
   12a59:	c9                   	leave  
   12a5a:	c3                   	ret    

0000000000012a5b <_Z7paddingPP10particle_tPii>:

void padding(particle_t** particles, int* nParticles, int minimalGrid) {
   12a5b:	f3 0f 1e fa          	endbr64 
   12a5f:	55                   	push   %rbp
   12a60:	48 89 e5             	mov    %rsp,%rbp
   12a63:	48 81 ec e0 00 00 00 	sub    $0xe0,%rsp
   12a6a:	48 89 bd 38 ff ff ff 	mov    %rdi,-0xc8(%rbp)
   12a71:	48 89 b5 30 ff ff ff 	mov    %rsi,-0xd0(%rbp)
   12a78:	89 95 2c ff ff ff    	mov    %edx,-0xd4(%rbp)
   12a7e:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   12a85:	00 00 
   12a87:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   12a8b:	31 c0                	xor    %eax,%eax
    int oriSize = *nParticles;
   12a8d:	48 8b 85 30 ff ff ff 	mov    -0xd0(%rbp),%rax
   12a94:	8b 00                	mov    (%rax),%eax
   12a96:	89 85 40 ff ff ff    	mov    %eax,-0xc0(%rbp)
    int newSize = sizeAfterPadding(oriSize, minimalGrid);
   12a9c:	8b 95 2c ff ff ff    	mov    -0xd4(%rbp),%edx
   12aa2:	8b 85 40 ff ff ff    	mov    -0xc0(%rbp),%eax
   12aa8:	89 d6                	mov    %edx,%esi
   12aaa:	89 c7                	mov    %eax,%edi
   12aac:	e8 48 31 00 00       	call   15bf9 <_Z16sizeAfterPaddingii>
   12ab1:	89 85 44 ff ff ff    	mov    %eax,-0xbc(%rbp)
    particle_t *oriP, *newP;
    double *oriPos, *oriVel, *oriAcc, *oriFeat;
    double *newPos, *newVel, *newAcc, *newFeat;
    oriP = *particles;
   12ab7:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
   12abe:	48 8b 00             	mov    (%rax),%rax
   12ac1:	48 89 85 50 ff ff ff 	mov    %rax,-0xb0(%rbp)
    oriPos = oriP->position;
   12ac8:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   12acf:	48 8b 40 08          	mov    0x8(%rax),%rax
   12ad3:	48 89 85 58 ff ff ff 	mov    %rax,-0xa8(%rbp)
    oriVel = oriP->velocity;
   12ada:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   12ae1:	48 8b 40 10          	mov    0x10(%rax),%rax
   12ae5:	48 89 85 60 ff ff ff 	mov    %rax,-0xa0(%rbp)
    oriAcc = oriP->acceleration;
   12aec:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   12af3:	48 8b 40 18          	mov    0x18(%rax),%rax
   12af7:	48 89 85 68 ff ff ff 	mov    %rax,-0x98(%rbp)
    oriFeat = oriP->features;
   12afe:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   12b05:	48 8b 40 20          	mov    0x20(%rax),%rax
   12b09:	48 89 85 70 ff ff ff 	mov    %rax,-0x90(%rbp)

    int ndim = oriP->ndim;
   12b10:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   12b17:	8b 40 30             	mov    0x30(%rax),%eax
   12b1a:	89 85 48 ff ff ff    	mov    %eax,-0xb8(%rbp)
    int nfeat = oriP->nfeat;
   12b20:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   12b27:	8b 40 34             	mov    0x34(%rax),%eax
   12b2a:	89 85 4c ff ff ff    	mov    %eax,-0xb4(%rbp)

    newP = static_cast<particle_t*>(realloc(oriP, newSize*sizeof(particle_t)));
   12b30:	8b 85 44 ff ff ff    	mov    -0xbc(%rbp),%eax
   12b36:	48 98                	cltq   
   12b38:	48 c1 e0 06          	shl    $0x6,%rax
   12b3c:	48 89 c2             	mov    %rax,%rdx
   12b3f:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   12b46:	48 89 d6             	mov    %rdx,%rsi
   12b49:	48 89 c7             	mov    %rax,%rdi
   12b4c:	e8 7f 7d ff ff       	call   a8d0 <realloc@plt>
   12b51:	48 89 85 78 ff ff ff 	mov    %rax,-0x88(%rbp)
    newPos = static_cast<double*>(realloc(oriPos, newSize*ndim*sizeof(double)));
   12b58:	8b 85 44 ff ff ff    	mov    -0xbc(%rbp),%eax
   12b5e:	0f af 85 48 ff ff ff 	imul   -0xb8(%rbp),%eax
   12b65:	48 98                	cltq   
   12b67:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12b6e:	00 
   12b6f:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   12b76:	48 89 d6             	mov    %rdx,%rsi
   12b79:	48 89 c7             	mov    %rax,%rdi
   12b7c:	e8 4f 7d ff ff       	call   a8d0 <realloc@plt>
   12b81:	48 89 45 80          	mov    %rax,-0x80(%rbp)
    newVel = static_cast<double*>(realloc(oriVel, newSize*ndim*sizeof(double)));
   12b85:	8b 85 44 ff ff ff    	mov    -0xbc(%rbp),%eax
   12b8b:	0f af 85 48 ff ff ff 	imul   -0xb8(%rbp),%eax
   12b92:	48 98                	cltq   
   12b94:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12b9b:	00 
   12b9c:	48 8b 85 60 ff ff ff 	mov    -0xa0(%rbp),%rax
   12ba3:	48 89 d6             	mov    %rdx,%rsi
   12ba6:	48 89 c7             	mov    %rax,%rdi
   12ba9:	e8 22 7d ff ff       	call   a8d0 <realloc@plt>
   12bae:	48 89 45 88          	mov    %rax,-0x78(%rbp)
    newAcc = static_cast<double*>(realloc(oriAcc, newSize*ndim*sizeof(double)));
   12bb2:	8b 85 44 ff ff ff    	mov    -0xbc(%rbp),%eax
   12bb8:	0f af 85 48 ff ff ff 	imul   -0xb8(%rbp),%eax
   12bbf:	48 98                	cltq   
   12bc1:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12bc8:	00 
   12bc9:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   12bd0:	48 89 d6             	mov    %rdx,%rsi
   12bd3:	48 89 c7             	mov    %rax,%rdi
   12bd6:	e8 f5 7c ff ff       	call   a8d0 <realloc@plt>
   12bdb:	48 89 45 90          	mov    %rax,-0x70(%rbp)
    newFeat = static_cast<double*>(realloc(oriFeat, newSize*nfeat*sizeof(double)));
   12bdf:	8b 85 44 ff ff ff    	mov    -0xbc(%rbp),%eax
   12be5:	0f af 85 4c ff ff ff 	imul   -0xb4(%rbp),%eax
   12bec:	48 98                	cltq   
   12bee:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12bf5:	00 
   12bf6:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   12bfd:	48 89 d6             	mov    %rdx,%rsi
   12c00:	48 89 c7             	mov    %rax,%rdi
   12c03:	e8 c8 7c ff ff       	call   a8d0 <realloc@plt>
   12c08:	48 89 45 98          	mov    %rax,-0x68(%rbp)

    assert(newP != nullptr);    assert(newPos != nullptr);      
   12c0c:	48 83 bd 78 ff ff ff 	cmpq   $0x0,-0x88(%rbp)
   12c13:	00 
   12c14:	75 28                	jne    12c3e <_Z7paddingPP10particle_tPii+0x1e3>
   12c16:	48 8d 05 33 41 00 00 	lea    0x4133(%rip),%rax        # 16d50 <_ZN3MPIL11LOCK_SHAREDE+0xb0>
   12c1d:	48 89 c1             	mov    %rax,%rcx
   12c20:	ba 3c 00 00 00       	mov    $0x3c,%edx
   12c25:	48 8d 05 ac 40 00 00 	lea    0x40ac(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12c2c:	48 89 c6             	mov    %rax,%rsi
   12c2f:	48 8d 05 40 41 00 00 	lea    0x4140(%rip),%rax        # 16d76 <_ZN3MPIL11LOCK_SHAREDE+0xd6>
   12c36:	48 89 c7             	mov    %rax,%rdi
   12c39:	e8 12 75 ff ff       	call   a150 <__assert_fail@plt>
   12c3e:	48 83 7d 80 00       	cmpq   $0x0,-0x80(%rbp)
   12c43:	75 28                	jne    12c6d <_Z7paddingPP10particle_tPii+0x212>
   12c45:	48 8d 05 04 41 00 00 	lea    0x4104(%rip),%rax        # 16d50 <_ZN3MPIL11LOCK_SHAREDE+0xb0>
   12c4c:	48 89 c1             	mov    %rax,%rcx
   12c4f:	ba 3c 00 00 00       	mov    $0x3c,%edx
   12c54:	48 8d 05 7d 40 00 00 	lea    0x407d(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12c5b:	48 89 c6             	mov    %rax,%rsi
   12c5e:	48 8d 05 21 41 00 00 	lea    0x4121(%rip),%rax        # 16d86 <_ZN3MPIL11LOCK_SHAREDE+0xe6>
   12c65:	48 89 c7             	mov    %rax,%rdi
   12c68:	e8 e3 74 ff ff       	call   a150 <__assert_fail@plt>
    assert(newVel != nullptr);  assert(newAcc != nullptr);
   12c6d:	48 83 7d 88 00       	cmpq   $0x0,-0x78(%rbp)
   12c72:	75 28                	jne    12c9c <_Z7paddingPP10particle_tPii+0x241>
   12c74:	48 8d 05 d5 40 00 00 	lea    0x40d5(%rip),%rax        # 16d50 <_ZN3MPIL11LOCK_SHAREDE+0xb0>
   12c7b:	48 89 c1             	mov    %rax,%rcx
   12c7e:	ba 3d 00 00 00       	mov    $0x3d,%edx
   12c83:	48 8d 05 4e 40 00 00 	lea    0x404e(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12c8a:	48 89 c6             	mov    %rax,%rsi
   12c8d:	48 8d 05 04 41 00 00 	lea    0x4104(%rip),%rax        # 16d98 <_ZN3MPIL11LOCK_SHAREDE+0xf8>
   12c94:	48 89 c7             	mov    %rax,%rdi
   12c97:	e8 b4 74 ff ff       	call   a150 <__assert_fail@plt>
   12c9c:	48 83 7d 90 00       	cmpq   $0x0,-0x70(%rbp)
   12ca1:	75 28                	jne    12ccb <_Z7paddingPP10particle_tPii+0x270>
   12ca3:	48 8d 05 a6 40 00 00 	lea    0x40a6(%rip),%rax        # 16d50 <_ZN3MPIL11LOCK_SHAREDE+0xb0>
   12caa:	48 89 c1             	mov    %rax,%rcx
   12cad:	ba 3d 00 00 00       	mov    $0x3d,%edx
   12cb2:	48 8d 05 1f 40 00 00 	lea    0x401f(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12cb9:	48 89 c6             	mov    %rax,%rsi
   12cbc:	48 8d 05 e7 40 00 00 	lea    0x40e7(%rip),%rax        # 16daa <_ZN3MPIL11LOCK_SHAREDE+0x10a>
   12cc3:	48 89 c7             	mov    %rax,%rdi
   12cc6:	e8 85 74 ff ff       	call   a150 <__assert_fail@plt>

    #pragma omp parallel for collapse(2)
   12ccb:	8b 85 48 ff ff ff    	mov    -0xb8(%rbp),%eax
   12cd1:	89 45 c0             	mov    %eax,-0x40(%rbp)
   12cd4:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   12cd8:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
   12cdc:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   12ce0:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
   12ce4:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   12ce8:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
   12cec:	8b 85 44 ff ff ff    	mov    -0xbc(%rbp),%eax
   12cf2:	89 45 bc             	mov    %eax,-0x44(%rbp)
   12cf5:	8b 85 40 ff ff ff    	mov    -0xc0(%rbp),%eax
   12cfb:	89 45 b8             	mov    %eax,-0x48(%rbp)
   12cfe:	48 8d 45 a0          	lea    -0x60(%rbp),%rax
   12d02:	b9 00 00 00 00       	mov    $0x0,%ecx
   12d07:	ba 00 00 00 00       	mov    $0x0,%edx
   12d0c:	48 89 c6             	mov    %rax,%rsi
   12d0f:	48 8d 05 f5 1e 00 00 	lea    0x1ef5(%rip),%rax        # 14c0b <_Z7paddingPP10particle_tPii._omp_fn.0>
   12d16:	48 89 c7             	mov    %rax,%rdi
   12d19:	e8 92 7a ff ff       	call   a7b0 <GOMP_parallel@plt>
            *(newPos+n*ndim+nn) = 0;
            *(newVel+n*ndim+nn) = 0;
            *(newAcc+n*ndim+nn) = 0;
        }
    
    #pragma omp parallel for collapse(2)
   12d1e:	8b 85 4c ff ff ff    	mov    -0xb4(%rbp),%eax
   12d24:	89 45 b0             	mov    %eax,-0x50(%rbp)
   12d27:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   12d2b:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
   12d2f:	8b 85 44 ff ff ff    	mov    -0xbc(%rbp),%eax
   12d35:	89 45 ac             	mov    %eax,-0x54(%rbp)
   12d38:	8b 85 40 ff ff ff    	mov    -0xc0(%rbp),%eax
   12d3e:	89 45 a8             	mov    %eax,-0x58(%rbp)
   12d41:	48 8d 45 a0          	lea    -0x60(%rbp),%rax
   12d45:	b9 00 00 00 00       	mov    $0x0,%ecx
   12d4a:	ba 00 00 00 00       	mov    $0x0,%edx
   12d4f:	48 89 c6             	mov    %rax,%rsi
   12d52:	48 8d 05 52 20 00 00 	lea    0x2052(%rip),%rax        # 14dab <_Z7paddingPP10particle_tPii._omp_fn.1>
   12d59:	48 89 c7             	mov    %rax,%rdi
   12d5c:	e8 4f 7a ff ff       	call   a7b0 <GOMP_parallel@plt>
    for(int n=oriSize; n<newSize; n++)
        for(int nn=0; nn<nfeat; nn++)
            *(newFeat+n*nfeat+nn) = 0;
    
    #pragma omp parallel for
   12d61:	8b 85 4c ff ff ff    	mov    -0xb4(%rbp),%eax
   12d67:	89 45 f4             	mov    %eax,-0xc(%rbp)
   12d6a:	8b 85 48 ff ff ff    	mov    -0xb8(%rbp),%eax
   12d70:	89 45 f0             	mov    %eax,-0x10(%rbp)
   12d73:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   12d77:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   12d7b:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   12d7f:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   12d83:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   12d87:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   12d8b:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   12d8f:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   12d93:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   12d9a:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   12d9e:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   12da5:	48 89 45 b8          	mov    %rax,-0x48(%rbp)
   12da9:	48 8b 85 60 ff ff ff 	mov    -0xa0(%rbp),%rax
   12db0:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
   12db4:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   12dbb:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
   12dbf:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   12dc6:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
   12dca:	8b 85 44 ff ff ff    	mov    -0xbc(%rbp),%eax
   12dd0:	89 45 ec             	mov    %eax,-0x14(%rbp)
   12dd3:	8b 85 40 ff ff ff    	mov    -0xc0(%rbp),%eax
   12dd9:	89 45 e8             	mov    %eax,-0x18(%rbp)
   12ddc:	48 8d 45 a0          	lea    -0x60(%rbp),%rax
   12de0:	b9 00 00 00 00       	mov    $0x0,%ecx
   12de5:	ba 00 00 00 00       	mov    $0x0,%edx
   12dea:	48 89 c6             	mov    %rax,%rsi
   12ded:	48 8d 05 db 20 00 00 	lea    0x20db(%rip),%rax        # 14ecf <_Z7paddingPP10particle_tPii._omp_fn.2>
   12df4:	48 89 c7             	mov    %rax,%rdi
   12df7:	e8 b4 79 ff ff       	call   a7b0 <GOMP_parallel@plt>
            (newP+n)->nfeat = nfeat;
            (newP+n)->updateAcceleration = newP->updateAcceleration;
        }
    }

    *particles = newP;
   12dfc:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
   12e03:	48 8b 95 78 ff ff ff 	mov    -0x88(%rbp),%rdx
   12e0a:	48 89 10             	mov    %rdx,(%rax)
    *nParticles = newSize;
   12e0d:	48 8b 85 30 ff ff ff 	mov    -0xd0(%rbp),%rax
   12e14:	8b 95 44 ff ff ff    	mov    -0xbc(%rbp),%edx
   12e1a:	89 10                	mov    %edx,(%rax)
}
   12e1c:	90                   	nop
   12e1d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12e21:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   12e28:	00 00 
   12e2a:	74 05                	je     12e31 <_Z7paddingPP10particle_tPii+0x3d6>
   12e2c:	e8 af 76 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   12e31:	c9                   	leave  
   12e32:	c3                   	ret    

0000000000012e33 <_Z7paddingP17chunk_particles_ti>:

void padding(chunk_particles_t* chunkParticles, int minimalGrid) {
   12e33:	f3 0f 1e fa          	endbr64 
   12e37:	55                   	push   %rbp
   12e38:	48 89 e5             	mov    %rsp,%rbp
   12e3b:	48 83 ec 30          	sub    $0x30,%rsp
   12e3f:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   12e43:	89 75 d4             	mov    %esi,-0x2c(%rbp)
   12e46:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   12e4d:	00 00 
   12e4f:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   12e53:	31 c0                	xor    %eax,%eax
    particle_t* particles = chunkParticles->particles;
   12e55:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12e59:	48 8b 00             	mov    (%rax),%rax
   12e5c:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
    int nParticles = chunkParticles->nParticle;
   12e60:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12e64:	8b 40 08             	mov    0x8(%rax),%eax
   12e67:	89 45 ec             	mov    %eax,-0x14(%rbp)
    padding(&particles, &nParticles, minimalGrid);
   12e6a:	8b 55 d4             	mov    -0x2c(%rbp),%edx
   12e6d:	48 8d 4d ec          	lea    -0x14(%rbp),%rcx
   12e71:	48 8d 45 f0          	lea    -0x10(%rbp),%rax
   12e75:	48 89 ce             	mov    %rcx,%rsi
   12e78:	48 89 c7             	mov    %rax,%rdi
   12e7b:	e8 db fb ff ff       	call   12a5b <_Z7paddingPP10particle_tPii>
    chunkParticles->particles = particles;
   12e80:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   12e84:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12e88:	48 89 10             	mov    %rdx,(%rax)
    chunkParticles->nParticle = nParticles;
   12e8b:	8b 55 ec             	mov    -0x14(%rbp),%edx
   12e8e:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12e92:	89 50 08             	mov    %edx,0x8(%rax)
}
   12e95:	90                   	nop
   12e96:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12e9a:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   12ea1:	00 00 
   12ea3:	74 05                	je     12eaa <_Z7paddingP17chunk_particles_ti+0x77>
   12ea5:	e8 36 76 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   12eaa:	c9                   	leave  
   12eab:	c3                   	ret    

0000000000012eac <_Z6shrinkPP10particle_tPi>:

void shrink(particle_t** particles, int* nParticles) {
   12eac:	f3 0f 1e fa          	endbr64 
   12eb0:	55                   	push   %rbp
   12eb1:	48 89 e5             	mov    %rsp,%rbp
   12eb4:	48 81 ec 90 00 00 00 	sub    $0x90,%rsp
   12ebb:	48 89 bd 78 ff ff ff 	mov    %rdi,-0x88(%rbp)
   12ec2:	48 89 b5 70 ff ff ff 	mov    %rsi,-0x90(%rbp)
    particle_t* oriParticles = *particles;
   12ec9:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   12ed0:	48 8b 00             	mov    (%rax),%rax
   12ed3:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
    int oriNParticles = *nParticles;
   12ed7:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   12ede:	8b 00                	mov    (%rax),%eax
   12ee0:	89 45 94             	mov    %eax,-0x6c(%rbp)
    int newSize;

    particle_t* ppt = oriParticles;
   12ee3:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12ee7:	48 89 45 98          	mov    %rax,-0x68(%rbp)
    double* dpt;
    for(newSize=0; newSize<oriNParticles; newSize++) {
   12eeb:	c7 45 80 00 00 00 00 	movl   $0x0,-0x80(%rbp)
   12ef2:	eb 18                	jmp    12f0c <_Z6shrinkPP10particle_tPi+0x60>
        if (!ppt->enabled) {
   12ef4:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   12ef8:	0f b6 40 38          	movzbl 0x38(%rax),%eax
   12efc:	83 f0 01             	xor    $0x1,%eax
   12eff:	84 c0                	test   %al,%al
   12f01:	75 13                	jne    12f16 <_Z6shrinkPP10particle_tPi+0x6a>
            break;
        }
        ppt++;
   12f03:	48 83 45 98 40       	addq   $0x40,-0x68(%rbp)
    for(newSize=0; newSize<oriNParticles; newSize++) {
   12f08:	83 45 80 01          	addl   $0x1,-0x80(%rbp)
   12f0c:	8b 45 80             	mov    -0x80(%rbp),%eax
   12f0f:	3b 45 94             	cmp    -0x6c(%rbp),%eax
   12f12:	7c e0                	jl     12ef4 <_Z6shrinkPP10particle_tPi+0x48>
   12f14:	eb 01                	jmp    12f17 <_Z6shrinkPP10particle_tPi+0x6b>
            break;
   12f16:	90                   	nop
    }
    double *oriPos, *oriVel, *oriAcc, *oriFeat;
    double *newPos, *newVel, *newAcc, *newFeat; 
    particle_t *oriP, *newP;

    oriPos = oriParticles->position;
   12f17:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12f1b:	48 8b 40 08          	mov    0x8(%rax),%rax
   12f1f:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
    newPos = static_cast<double*>(realloc(oriPos, newSize*oriParticles->ndim*sizeof(double)));
   12f23:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12f27:	8b 40 30             	mov    0x30(%rax),%eax
   12f2a:	0f af 45 80          	imul   -0x80(%rbp),%eax
   12f2e:	48 98                	cltq   
   12f30:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12f37:	00 
   12f38:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   12f3c:	48 89 d6             	mov    %rdx,%rsi
   12f3f:	48 89 c7             	mov    %rax,%rdi
   12f42:	e8 89 79 ff ff       	call   a8d0 <realloc@plt>
   12f47:	48 89 45 b8          	mov    %rax,-0x48(%rbp)
    assert(newPos != nullptr);
   12f4b:	48 83 7d b8 00       	cmpq   $0x0,-0x48(%rbp)
   12f50:	75 28                	jne    12f7a <_Z6shrinkPP10particle_tPi+0xce>
   12f52:	48 8d 05 67 3e 00 00 	lea    0x3e67(%rip),%rax        # 16dc0 <_ZN3MPIL11LOCK_SHAREDE+0x120>
   12f59:	48 89 c1             	mov    %rax,%rcx
   12f5c:	ba 7e 00 00 00       	mov    $0x7e,%edx
   12f61:	48 8d 05 70 3d 00 00 	lea    0x3d70(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12f68:	48 89 c6             	mov    %rax,%rsi
   12f6b:	48 8d 05 14 3e 00 00 	lea    0x3e14(%rip),%rax        # 16d86 <_ZN3MPIL11LOCK_SHAREDE+0xe6>
   12f72:	48 89 c7             	mov    %rax,%rdi
   12f75:	e8 d6 71 ff ff       	call   a150 <__assert_fail@plt>

    oriVel = oriParticles->velocity;
   12f7a:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12f7e:	48 8b 40 10          	mov    0x10(%rax),%rax
   12f82:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
    newVel = static_cast<double*>(realloc(oriVel, newSize*oriParticles->ndim*sizeof(double)));
   12f86:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12f8a:	8b 40 30             	mov    0x30(%rax),%eax
   12f8d:	0f af 45 80          	imul   -0x80(%rbp),%eax
   12f91:	48 98                	cltq   
   12f93:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12f9a:	00 
   12f9b:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   12f9f:	48 89 d6             	mov    %rdx,%rsi
   12fa2:	48 89 c7             	mov    %rax,%rdi
   12fa5:	e8 26 79 ff ff       	call   a8d0 <realloc@plt>
   12faa:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
    assert(newVel != nullptr);
   12fae:	48 83 7d c8 00       	cmpq   $0x0,-0x38(%rbp)
   12fb3:	75 28                	jne    12fdd <_Z6shrinkPP10particle_tPi+0x131>
   12fb5:	48 8d 05 04 3e 00 00 	lea    0x3e04(%rip),%rax        # 16dc0 <_ZN3MPIL11LOCK_SHAREDE+0x120>
   12fbc:	48 89 c1             	mov    %rax,%rcx
   12fbf:	ba 82 00 00 00       	mov    $0x82,%edx
   12fc4:	48 8d 05 0d 3d 00 00 	lea    0x3d0d(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12fcb:	48 89 c6             	mov    %rax,%rsi
   12fce:	48 8d 05 c3 3d 00 00 	lea    0x3dc3(%rip),%rax        # 16d98 <_ZN3MPIL11LOCK_SHAREDE+0xf8>
   12fd5:	48 89 c7             	mov    %rax,%rdi
   12fd8:	e8 73 71 ff ff       	call   a150 <__assert_fail@plt>
    
    oriAcc = oriParticles->acceleration;
   12fdd:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12fe1:	48 8b 40 18          	mov    0x18(%rax),%rax
   12fe5:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
    newAcc = static_cast<double*>(realloc(oriAcc, newSize*oriParticles->ndim*sizeof(double)));
   12fe9:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12fed:	8b 40 30             	mov    0x30(%rax),%eax
   12ff0:	0f af 45 80          	imul   -0x80(%rbp),%eax
   12ff4:	48 98                	cltq   
   12ff6:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12ffd:	00 
   12ffe:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   13002:	48 89 d6             	mov    %rdx,%rsi
   13005:	48 89 c7             	mov    %rax,%rdi
   13008:	e8 c3 78 ff ff       	call   a8d0 <realloc@plt>
   1300d:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
    assert(newAcc != nullptr);
   13011:	48 83 7d d8 00       	cmpq   $0x0,-0x28(%rbp)
   13016:	75 28                	jne    13040 <_Z6shrinkPP10particle_tPi+0x194>
   13018:	48 8d 05 a1 3d 00 00 	lea    0x3da1(%rip),%rax        # 16dc0 <_ZN3MPIL11LOCK_SHAREDE+0x120>
   1301f:	48 89 c1             	mov    %rax,%rcx
   13022:	ba 86 00 00 00       	mov    $0x86,%edx
   13027:	48 8d 05 aa 3c 00 00 	lea    0x3caa(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   1302e:	48 89 c6             	mov    %rax,%rsi
   13031:	48 8d 05 72 3d 00 00 	lea    0x3d72(%rip),%rax        # 16daa <_ZN3MPIL11LOCK_SHAREDE+0x10a>
   13038:	48 89 c7             	mov    %rax,%rdi
   1303b:	e8 10 71 ff ff       	call   a150 <__assert_fail@plt>

    oriFeat = oriParticles->features;
   13040:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   13044:	48 8b 40 20          	mov    0x20(%rax),%rax
   13048:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
    newFeat = static_cast<double*>(realloc(oriFeat, newSize*oriParticles->nfeat*sizeof(double)));
   1304c:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   13050:	8b 40 34             	mov    0x34(%rax),%eax
   13053:	0f af 45 80          	imul   -0x80(%rbp),%eax
   13057:	48 98                	cltq   
   13059:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   13060:	00 
   13061:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   13065:	48 89 d6             	mov    %rdx,%rsi
   13068:	48 89 c7             	mov    %rax,%rdi
   1306b:	e8 60 78 ff ff       	call   a8d0 <realloc@plt>
   13070:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    assert(newFeat != nullptr);
   13074:	48 83 7d e8 00       	cmpq   $0x0,-0x18(%rbp)
   13079:	75 28                	jne    130a3 <_Z6shrinkPP10particle_tPi+0x1f7>
   1307b:	48 8d 05 3e 3d 00 00 	lea    0x3d3e(%rip),%rax        # 16dc0 <_ZN3MPIL11LOCK_SHAREDE+0x120>
   13082:	48 89 c1             	mov    %rax,%rcx
   13085:	ba 8a 00 00 00       	mov    $0x8a,%edx
   1308a:	48 8d 05 47 3c 00 00 	lea    0x3c47(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   13091:	48 89 c6             	mov    %rax,%rsi
   13094:	48 8d 05 45 3d 00 00 	lea    0x3d45(%rip),%rax        # 16de0 <_ZN3MPIL11LOCK_SHAREDE+0x140>
   1309b:	48 89 c7             	mov    %rax,%rdi
   1309e:	e8 ad 70 ff ff       	call   a150 <__assert_fail@plt>

    oriP = oriParticles;
   130a3:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   130a7:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
    newP = static_cast<particle_t*>(realloc(oriP, newSize*sizeof(particle_t)));
   130ab:	8b 45 80             	mov    -0x80(%rbp),%eax
   130ae:	48 98                	cltq   
   130b0:	48 c1 e0 06          	shl    $0x6,%rax
   130b4:	48 89 c2             	mov    %rax,%rdx
   130b7:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   130bb:	48 89 d6             	mov    %rdx,%rsi
   130be:	48 89 c7             	mov    %rax,%rdi
   130c1:	e8 0a 78 ff ff       	call   a8d0 <realloc@plt>
   130c6:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    assert(newP != nullptr);
   130ca:	48 83 7d f8 00       	cmpq   $0x0,-0x8(%rbp)
   130cf:	75 28                	jne    130f9 <_Z6shrinkPP10particle_tPi+0x24d>
   130d1:	48 8d 05 e8 3c 00 00 	lea    0x3ce8(%rip),%rax        # 16dc0 <_ZN3MPIL11LOCK_SHAREDE+0x120>
   130d8:	48 89 c1             	mov    %rax,%rcx
   130db:	ba 8e 00 00 00       	mov    $0x8e,%edx
   130e0:	48 8d 05 f1 3b 00 00 	lea    0x3bf1(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   130e7:	48 89 c6             	mov    %rax,%rsi
   130ea:	48 8d 05 85 3c 00 00 	lea    0x3c85(%rip),%rax        # 16d76 <_ZN3MPIL11LOCK_SHAREDE+0xd6>
   130f1:	48 89 c7             	mov    %rax,%rdi
   130f4:	e8 57 70 ff ff       	call   a150 <__assert_fail@plt>

    if(oriPos != newPos) {
   130f9:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   130fd:	48 3b 45 b8          	cmp    -0x48(%rbp),%rax
   13101:	74 47                	je     1314a <_Z6shrinkPP10particle_tPi+0x29e>
        ppt = newP;
   13103:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   13107:	48 89 45 98          	mov    %rax,-0x68(%rbp)
        dpt = newPos;
   1310b:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   1310f:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
        for(int n=0; n<newSize; n++) {
   13113:	c7 45 84 00 00 00 00 	movl   $0x0,-0x7c(%rbp)
   1311a:	eb 26                	jmp    13142 <_Z6shrinkPP10particle_tPi+0x296>
            ppt->position = dpt;
   1311c:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   13120:	48 8b 55 a0          	mov    -0x60(%rbp),%rdx
   13124:	48 89 50 08          	mov    %rdx,0x8(%rax)
            dpt+=ppt->ndim;
   13128:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1312c:	8b 40 30             	mov    0x30(%rax),%eax
   1312f:	48 98                	cltq   
   13131:	48 c1 e0 03          	shl    $0x3,%rax
   13135:	48 01 45 a0          	add    %rax,-0x60(%rbp)
            ppt++;
   13139:	48 83 45 98 40       	addq   $0x40,-0x68(%rbp)
        for(int n=0; n<newSize; n++) {
   1313e:	83 45 84 01          	addl   $0x1,-0x7c(%rbp)
   13142:	8b 45 84             	mov    -0x7c(%rbp),%eax
   13145:	3b 45 80             	cmp    -0x80(%rbp),%eax
   13148:	7c d2                	jl     1311c <_Z6shrinkPP10particle_tPi+0x270>
        }
    }

    if(oriVel != newVel) {
   1314a:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   1314e:	48 3b 45 c8          	cmp    -0x38(%rbp),%rax
   13152:	74 47                	je     1319b <_Z6shrinkPP10particle_tPi+0x2ef>
        ppt = newP;
   13154:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   13158:	48 89 45 98          	mov    %rax,-0x68(%rbp)
        dpt = newVel;
   1315c:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   13160:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
        for(int n=0; n<newSize; n++) {
   13164:	c7 45 88 00 00 00 00 	movl   $0x0,-0x78(%rbp)
   1316b:	eb 26                	jmp    13193 <_Z6shrinkPP10particle_tPi+0x2e7>
            ppt->velocity = dpt;
   1316d:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   13171:	48 8b 55 a0          	mov    -0x60(%rbp),%rdx
   13175:	48 89 50 10          	mov    %rdx,0x10(%rax)
            dpt+=ppt->ndim;
   13179:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1317d:	8b 40 30             	mov    0x30(%rax),%eax
   13180:	48 98                	cltq   
   13182:	48 c1 e0 03          	shl    $0x3,%rax
   13186:	48 01 45 a0          	add    %rax,-0x60(%rbp)
            ppt++;
   1318a:	48 83 45 98 40       	addq   $0x40,-0x68(%rbp)
        for(int n=0; n<newSize; n++) {
   1318f:	83 45 88 01          	addl   $0x1,-0x78(%rbp)
   13193:	8b 45 88             	mov    -0x78(%rbp),%eax
   13196:	3b 45 80             	cmp    -0x80(%rbp),%eax
   13199:	7c d2                	jl     1316d <_Z6shrinkPP10particle_tPi+0x2c1>
        }
    }

    if(oriAcc != newAcc) {
   1319b:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   1319f:	48 3b 45 d8          	cmp    -0x28(%rbp),%rax
   131a3:	74 47                	je     131ec <_Z6shrinkPP10particle_tPi+0x340>
        ppt = newP;
   131a5:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   131a9:	48 89 45 98          	mov    %rax,-0x68(%rbp)
        dpt = newAcc;
   131ad:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   131b1:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
        for(int n=0; n<newSize; n++) {
   131b5:	c7 45 8c 00 00 00 00 	movl   $0x0,-0x74(%rbp)
   131bc:	eb 26                	jmp    131e4 <_Z6shrinkPP10particle_tPi+0x338>
            ppt->acceleration = dpt;
   131be:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   131c2:	48 8b 55 a0          	mov    -0x60(%rbp),%rdx
   131c6:	48 89 50 18          	mov    %rdx,0x18(%rax)
            dpt+=ppt->ndim;
   131ca:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   131ce:	8b 40 30             	mov    0x30(%rax),%eax
   131d1:	48 98                	cltq   
   131d3:	48 c1 e0 03          	shl    $0x3,%rax
   131d7:	48 01 45 a0          	add    %rax,-0x60(%rbp)
            ppt++;
   131db:	48 83 45 98 40       	addq   $0x40,-0x68(%rbp)
        for(int n=0; n<newSize; n++) {
   131e0:	83 45 8c 01          	addl   $0x1,-0x74(%rbp)
   131e4:	8b 45 8c             	mov    -0x74(%rbp),%eax
   131e7:	3b 45 80             	cmp    -0x80(%rbp),%eax
   131ea:	7c d2                	jl     131be <_Z6shrinkPP10particle_tPi+0x312>
        }
    }

    if(oriFeat != newFeat) {
   131ec:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   131f0:	48 3b 45 e8          	cmp    -0x18(%rbp),%rax
   131f4:	74 47                	je     1323d <_Z6shrinkPP10particle_tPi+0x391>
        ppt = newP;
   131f6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   131fa:	48 89 45 98          	mov    %rax,-0x68(%rbp)
        dpt = newFeat;
   131fe:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   13202:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
        for(int n=0; n<newSize; n++) {
   13206:	c7 45 90 00 00 00 00 	movl   $0x0,-0x70(%rbp)
   1320d:	eb 26                	jmp    13235 <_Z6shrinkPP10particle_tPi+0x389>
            ppt->features = dpt;
   1320f:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   13213:	48 8b 55 a0          	mov    -0x60(%rbp),%rdx
   13217:	48 89 50 20          	mov    %rdx,0x20(%rax)
            dpt+=ppt->nfeat;
   1321b:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1321f:	8b 40 34             	mov    0x34(%rax),%eax
   13222:	48 98                	cltq   
   13224:	48 c1 e0 03          	shl    $0x3,%rax
   13228:	48 01 45 a0          	add    %rax,-0x60(%rbp)
            ppt++;
   1322c:	48 83 45 98 40       	addq   $0x40,-0x68(%rbp)
        for(int n=0; n<newSize; n++) {
   13231:	83 45 90 01          	addl   $0x1,-0x70(%rbp)
   13235:	8b 45 90             	mov    -0x70(%rbp),%eax
   13238:	3b 45 80             	cmp    -0x80(%rbp),%eax
   1323b:	7c d2                	jl     1320f <_Z6shrinkPP10particle_tPi+0x363>
        }
    }

    *particles = newP;
   1323d:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13244:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   13248:	48 89 10             	mov    %rdx,(%rax)
    *nParticles = newSize;
   1324b:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   13252:	8b 55 80             	mov    -0x80(%rbp),%edx
   13255:	89 10                	mov    %edx,(%rax)
}
   13257:	90                   	nop
   13258:	c9                   	leave  
   13259:	c3                   	ret    

000000000001325a <_Z6shrinkP17chunk_particles_t>:

void shrink(chunk_particles_t* chunkParticles) {
   1325a:	f3 0f 1e fa          	endbr64 
   1325e:	55                   	push   %rbp
   1325f:	48 89 e5             	mov    %rsp,%rbp
   13262:	48 83 ec 30          	sub    $0x30,%rsp
   13266:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   1326a:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   13271:	00 00 
   13273:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   13277:	31 c0                	xor    %eax,%eax
    particle_t* particles = chunkParticles->particles;
   13279:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   1327d:	48 8b 00             	mov    (%rax),%rax
   13280:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
    int nParticle = chunkParticles->nParticle;
   13284:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   13288:	8b 40 08             	mov    0x8(%rax),%eax
   1328b:	89 45 ec             	mov    %eax,-0x14(%rbp)
    shrink(&particles, &nParticle);
   1328e:	48 8d 55 ec          	lea    -0x14(%rbp),%rdx
   13292:	48 8d 45 f0          	lea    -0x10(%rbp),%rax
   13296:	48 89 d6             	mov    %rdx,%rsi
   13299:	48 89 c7             	mov    %rax,%rdi
   1329c:	e8 0b fc ff ff       	call   12eac <_Z6shrinkPP10particle_tPi>
    chunkParticles->particles = particles;
   132a1:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   132a5:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   132a9:	48 89 10             	mov    %rdx,(%rax)
    chunkParticles->nParticle = nParticle;
   132ac:	8b 55 ec             	mov    -0x14(%rbp),%edx
   132af:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   132b3:	89 50 08             	mov    %edx,0x8(%rax)
}
   132b6:	90                   	nop
   132b7:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   132bb:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   132c2:	00 00 
   132c4:	74 05                	je     132cb <_Z6shrinkP17chunk_particles_t+0x71>
   132c6:	e8 15 72 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   132cb:	c9                   	leave  
   132cc:	c3                   	ret    

00000000000132cd <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t>:

void scatter(particle_t** localParticlesX, int* localNParticlesX,
            particle_t** localParticlesY, int* localNParticlesY,
            particle_t* wholeParticles, int wholeNParticles, 
            topology_t* topology) {
   132cd:	f3 0f 1e fa          	endbr64 
   132d1:	55                   	push   %rbp
   132d2:	48 89 e5             	mov    %rsp,%rbp
   132d5:	48 81 ec e0 00 00 00 	sub    $0xe0,%rsp
   132dc:	48 89 bd 58 ff ff ff 	mov    %rdi,-0xa8(%rbp)
   132e3:	48 89 b5 50 ff ff ff 	mov    %rsi,-0xb0(%rbp)
   132ea:	48 89 95 48 ff ff ff 	mov    %rdx,-0xb8(%rbp)
   132f1:	48 89 8d 40 ff ff ff 	mov    %rcx,-0xc0(%rbp)
   132f8:	4c 89 85 38 ff ff ff 	mov    %r8,-0xc8(%rbp)
   132ff:	44 89 8d 34 ff ff ff 	mov    %r9d,-0xcc(%rbp)
   13306:	48 8b 45 10          	mov    0x10(%rbp),%rax
   1330a:	48 89 85 28 ff ff ff 	mov    %rax,-0xd8(%rbp)
   13311:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   13318:	00 00 
   1331a:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   1331e:	31 c0                	xor    %eax,%eax
    bool wholeParticlesInitNull = false;
   13320:	c6 85 6f ff ff ff 00 	movb   $0x0,-0x91(%rbp)
    if(wholeParticles==nullptr) {
   13327:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
   1332e:	48 85 c0             	test   %rax,%rax
   13331:	75 25                	jne    13358 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x8b>
        wholeParticlesInitNull = true;
   13333:	c6 85 6f ff ff ff 01 	movb   $0x1,-0x91(%rbp)
        alloc_particles(&wholeParticles, 1, 1, 1);
   1333a:	48 8d 85 38 ff ff ff 	lea    -0xc8(%rbp),%rax
   13341:	b9 01 00 00 00       	mov    $0x1,%ecx
   13346:	ba 01 00 00 00       	mov    $0x1,%edx
   1334b:	be 01 00 00 00       	mov    $0x1,%esi
   13350:	48 89 c7             	mov    %rax,%rdi
   13353:	e8 f7 ea ff ff       	call   11e4f <_Z15alloc_particlesPP10particle_tiii>
    }

    const int& rankX = topology->rankX;
   13358:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   1335f:	48 89 85 78 ff ff ff 	mov    %rax,-0x88(%rbp)
    const int& rankY = topology->rankY;
   13366:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   1336d:	48 83 c0 04          	add    $0x4,%rax
   13371:	48 89 45 80          	mov    %rax,-0x80(%rbp)
    const int& gridX = topology->gridX;
   13375:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   1337c:	48 83 c0 08          	add    $0x8,%rax
   13380:	48 89 45 88          	mov    %rax,-0x78(%rbp)
    const int& gridY = topology->gridY;
   13384:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   1338b:	48 83 c0 0c          	add    $0xc,%rax
   1338f:	48 89 45 90          	mov    %rax,-0x70(%rbp)
    if(rankX==0 && rankY==0) {
   13393:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   1339a:	8b 00                	mov    (%rax),%eax
   1339c:	85 c0                	test   %eax,%eax
   1339e:	0f 85 84 00 00 00    	jne    13428 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x15b>
   133a4:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   133a8:	8b 00                	mov    (%rax),%eax
   133aa:	85 c0                	test   %eax,%eax
   133ac:	75 7a                	jne    13428 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x15b>
        assert(wholeNParticles%gridX == 0);
   133ae:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   133b2:	8b 08                	mov    (%rax),%ecx
   133b4:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   133ba:	99                   	cltd   
   133bb:	f7 f9                	idiv   %ecx
   133bd:	89 d0                	mov    %edx,%eax
   133bf:	85 c0                	test   %eax,%eax
   133c1:	74 28                	je     133eb <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x11e>
   133c3:	48 8d 05 2e 3a 00 00 	lea    0x3a2e(%rip),%rax        # 16df8 <_ZN3MPIL11LOCK_SHAREDE+0x158>
   133ca:	48 89 c1             	mov    %rax,%rcx
   133cd:	ba d3 00 00 00       	mov    $0xd3,%edx
   133d2:	48 8d 05 ff 38 00 00 	lea    0x38ff(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   133d9:	48 89 c6             	mov    %rax,%rsi
   133dc:	48 8d 05 69 3a 00 00 	lea    0x3a69(%rip),%rax        # 16e4c <_ZN3MPIL11LOCK_SHAREDE+0x1ac>
   133e3:	48 89 c7             	mov    %rax,%rdi
   133e6:	e8 65 6d ff ff       	call   a150 <__assert_fail@plt>
        assert(wholeNParticles%gridY == 0);
   133eb:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   133ef:	8b 08                	mov    (%rax),%ecx
   133f1:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   133f7:	99                   	cltd   
   133f8:	f7 f9                	idiv   %ecx
   133fa:	89 d0                	mov    %edx,%eax
   133fc:	85 c0                	test   %eax,%eax
   133fe:	74 28                	je     13428 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x15b>
   13400:	48 8d 05 f1 39 00 00 	lea    0x39f1(%rip),%rax        # 16df8 <_ZN3MPIL11LOCK_SHAREDE+0x158>
   13407:	48 89 c1             	mov    %rax,%rcx
   1340a:	ba d4 00 00 00       	mov    $0xd4,%edx
   1340f:	48 8d 05 c2 38 00 00 	lea    0x38c2(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   13416:	48 89 c6             	mov    %rax,%rsi
   13419:	48 8d 05 47 3a 00 00 	lea    0x3a47(%rip),%rax        # 16e67 <_ZN3MPIL11LOCK_SHAREDE+0x1c7>
   13420:	48 89 c7             	mov    %rax,%rdi
   13423:	e8 28 6d ff ff       	call   a150 <__assert_fail@plt>
    }

    int* localParticlesMetadata = static_cast<int*>(malloc(4*sizeof(int)));
   13428:	bf 10 00 00 00       	mov    $0x10,%edi
   1342d:	e8 fe 71 ff ff       	call   a630 <malloc@plt>
   13432:	48 89 45 98          	mov    %rax,-0x68(%rbp)
    assert(localParticlesMetadata != nullptr);
   13436:	48 83 7d 98 00       	cmpq   $0x0,-0x68(%rbp)
   1343b:	75 28                	jne    13465 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x198>
   1343d:	48 8d 05 b4 39 00 00 	lea    0x39b4(%rip),%rax        # 16df8 <_ZN3MPIL11LOCK_SHAREDE+0x158>
   13444:	48 89 c1             	mov    %rax,%rcx
   13447:	ba d8 00 00 00       	mov    $0xd8,%edx
   1344c:	48 8d 05 85 38 00 00 	lea    0x3885(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   13453:	48 89 c6             	mov    %rax,%rsi
   13456:	48 8d 05 2b 3a 00 00 	lea    0x3a2b(%rip),%rax        # 16e88 <_ZN3MPIL11LOCK_SHAREDE+0x1e8>
   1345d:	48 89 c7             	mov    %rax,%rdi
   13460:	e8 eb 6c ff ff       	call   a150 <__assert_fail@plt>

    if(rankX==0 && rankY==0) {
   13465:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   1346c:	8b 00                	mov    (%rax),%eax
   1346e:	85 c0                	test   %eax,%eax
   13470:	75 62                	jne    134d4 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x207>
   13472:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   13476:	8b 00                	mov    (%rax),%eax
   13478:	85 c0                	test   %eax,%eax
   1347a:	75 58                	jne    134d4 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x207>
        *(localParticlesMetadata+0) = wholeNParticles / gridX;
   1347c:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   13480:	8b 38                	mov    (%rax),%edi
   13482:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   13488:	99                   	cltd   
   13489:	f7 ff                	idiv   %edi
   1348b:	89 c2                	mov    %eax,%edx
   1348d:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   13491:	89 10                	mov    %edx,(%rax)
        *(localParticlesMetadata+1) = wholeNParticles / gridY;
   13493:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   13497:	8b 38                	mov    (%rax),%edi
   13499:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1349d:	48 8d 48 04          	lea    0x4(%rax),%rcx
   134a1:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   134a7:	99                   	cltd   
   134a8:	f7 ff                	idiv   %edi
   134aa:	89 01                	mov    %eax,(%rcx)
        *(localParticlesMetadata+2) = wholeParticles->ndim;
   134ac:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
   134b3:	48 8b 55 98          	mov    -0x68(%rbp),%rdx
   134b7:	48 83 c2 08          	add    $0x8,%rdx
   134bb:	8b 40 30             	mov    0x30(%rax),%eax
   134be:	89 02                	mov    %eax,(%rdx)
        *(localParticlesMetadata+3) = wholeParticles->nfeat;
   134c0:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
   134c7:	48 8b 55 98          	mov    -0x68(%rbp),%rdx
   134cb:	48 83 c2 0c          	add    $0xc,%rdx
   134cf:	8b 40 34             	mov    0x34(%rax),%eax
   134d2:	89 02                	mov    %eax,(%rdx)
    }

    MPI_Bcast(localParticlesMetadata, 4, MPI_INT, 0, MPI_COMM_WORLD);
   134d4:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   134d8:	4c 8d 05 41 a5 00 00 	lea    0xa541(%rip),%r8        # 1da20 <ompi_mpi_comm_world@@Base>
   134df:	b9 00 00 00 00       	mov    $0x0,%ecx
   134e4:	48 8d 15 35 ab 00 00 	lea    0xab35(%rip),%rdx        # 1e020 <ompi_mpi_int@@Base>
   134eb:	be 04 00 00 00       	mov    $0x4,%esi
   134f0:	48 89 c7             	mov    %rax,%rdi
   134f3:	e8 78 69 ff ff       	call   9e70 <MPI_Bcast@plt>
    MPI_Barrier(MPI_COMM_WORLD);
   134f8:	48 8d 05 21 a5 00 00 	lea    0xa521(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   134ff:	48 89 c7             	mov    %rax,%rdi
   13502:	e8 d9 6c ff ff       	call   a1e0 <MPI_Barrier@plt>

    int ndim, nfeat;
    *localNParticlesX = *(localParticlesMetadata+0);
   13507:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1350b:	8b 10                	mov    (%rax),%edx
   1350d:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   13514:	89 10                	mov    %edx,(%rax)
    *localNParticlesY = *(localParticlesMetadata+1);
   13516:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1351a:	8b 50 04             	mov    0x4(%rax),%edx
   1351d:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   13524:	89 10                	mov    %edx,(%rax)
    ndim = *(localParticlesMetadata+2);
   13526:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1352a:	8b 40 08             	mov    0x8(%rax),%eax
   1352d:	89 85 70 ff ff ff    	mov    %eax,-0x90(%rbp)
    nfeat = *(localParticlesMetadata+3);
   13533:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   13537:	8b 40 0c             	mov    0xc(%rax),%eax
   1353a:	89 85 74 ff ff ff    	mov    %eax,-0x8c(%rbp)
    alloc_particles(localParticlesX, *localNParticlesX, ndim, nfeat);
   13540:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   13547:	8b 30                	mov    (%rax),%esi
   13549:	8b 8d 74 ff ff ff    	mov    -0x8c(%rbp),%ecx
   1354f:	8b 95 70 ff ff ff    	mov    -0x90(%rbp),%edx
   13555:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   1355c:	48 89 c7             	mov    %rax,%rdi
   1355f:	e8 eb e8 ff ff       	call   11e4f <_Z15alloc_particlesPP10particle_tiii>
    alloc_particles(localParticlesY, *localNParticlesY, ndim, nfeat);
   13564:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   1356b:	8b 30                	mov    (%rax),%esi
   1356d:	8b 8d 74 ff ff ff    	mov    -0x8c(%rbp),%ecx
   13573:	8b 95 70 ff ff ff    	mov    -0x90(%rbp),%edx
   13579:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13580:	48 89 c7             	mov    %rax,%rdi
   13583:	e8 c7 e8 ff ff       	call   11e4f <_Z15alloc_particlesPP10particle_tiii>

    double *pos, *vel, *acc, *feat;

    if(rankY==0) {
   13588:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   1358c:	8b 00                	mov    (%rax),%eax
   1358e:	85 c0                	test   %eax,%eax
   13590:	0f 85 cc 04 00 00    	jne    13a62 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x795>
        pos = (*localParticlesX)->position;
   13596:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   1359d:	48 8b 00             	mov    (%rax),%rax
   135a0:	48 8b 40 08          	mov    0x8(%rax),%rax
   135a4:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
        vel = (*localParticlesX)->velocity;
   135a8:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   135af:	48 8b 00             	mov    (%rax),%rax
   135b2:	48 8b 40 10          	mov    0x10(%rax),%rax
   135b6:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
        acc = (*localParticlesX)->acceleration;
   135ba:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   135c1:	48 8b 00             	mov    (%rax),%rax
   135c4:	48 8b 40 18          	mov    0x18(%rax),%rax
   135c8:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
        feat = (*localParticlesX)->features;
   135cc:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   135d3:	48 8b 00             	mov    (%rax),%rax
   135d6:	48 8b 40 20          	mov    0x20(%rax),%rax
   135da:	48 89 45 b8          	mov    %rax,-0x48(%rbp)
        printf("%s:%d, %lx %lx %lx %lx\n", __FILE__, __LINE__, pos, vel, acc, feat);
   135de:	48 8b 4d b0          	mov    -0x50(%rbp),%rcx
   135e2:	48 8b 55 a8          	mov    -0x58(%rbp),%rdx
   135e6:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   135ea:	48 83 ec 08          	sub    $0x8,%rsp
   135ee:	ff 75 b8             	push   -0x48(%rbp)
   135f1:	49 89 c9             	mov    %rcx,%r9
   135f4:	49 89 d0             	mov    %rdx,%r8
   135f7:	48 89 c1             	mov    %rax,%rcx
   135fa:	ba f3 00 00 00       	mov    $0xf3,%edx
   135ff:	48 8d 05 d2 36 00 00 	lea    0x36d2(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   13606:	48 89 c6             	mov    %rax,%rsi
   13609:	48 8d 05 9a 38 00 00 	lea    0x389a(%rip),%rax        # 16eaa <_ZN3MPIL11LOCK_SHAREDE+0x20a>
   13610:	48 89 c7             	mov    %rax,%rdi
   13613:	b8 00 00 00 00       	mov    $0x0,%eax
   13618:	e8 c3 67 ff ff       	call   9de0 <printf@plt>
   1361d:	48 83 c4 10          	add    $0x10,%rsp

        MPI_Barrier(topology->scatterXComm);
   13621:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13628:	48 8b 40 20          	mov    0x20(%rax),%rax
   1362c:	48 89 c7             	mov    %rax,%rdi
   1362f:	e8 ac 6b ff ff       	call   a1e0 <MPI_Barrier@plt>
        printf("(%d, %d) %s %d %d@%d\n", topology->rankX, topology->rankY, __FILE__, __LINE__, *localNParticlesX, wholeNParticles);
   13634:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   1363b:	8b 30                	mov    (%rax),%esi
   1363d:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13644:	8b 50 04             	mov    0x4(%rax),%edx
   13647:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   1364e:	8b 00                	mov    (%rax),%eax
   13650:	48 83 ec 08          	sub    $0x8,%rsp
   13654:	8b 8d 34 ff ff ff    	mov    -0xcc(%rbp),%ecx
   1365a:	51                   	push   %rcx
   1365b:	41 89 f1             	mov    %esi,%r9d
   1365e:	41 b8 f6 00 00 00    	mov    $0xf6,%r8d
   13664:	48 8d 0d 6d 36 00 00 	lea    0x366d(%rip),%rcx        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   1366b:	89 c6                	mov    %eax,%esi
   1366d:	48 8d 05 4e 38 00 00 	lea    0x384e(%rip),%rax        # 16ec2 <_ZN3MPIL11LOCK_SHAREDE+0x222>
   13674:	48 89 c7             	mov    %rax,%rdi
   13677:	b8 00 00 00 00       	mov    $0x0,%eax
   1367c:	e8 5f 67 ff ff       	call   9de0 <printf@plt>
   13681:	48 83 c4 10          	add    $0x10,%rsp
        MPI_Scatter(
   13685:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   1368c:	48 8b 70 20          	mov    0x20(%rax),%rsi
            wholeParticles,
            wholeNParticles*sizeof(particle_t)/gridX,
            MPI_BYTE,
            (*localParticlesX),
            (*localNParticlesX)*sizeof(particle_t),
   13690:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   13697:	8b 00                	mov    (%rax),%eax
   13699:	48 98                	cltq   
        MPI_Scatter(
   1369b:	c1 e0 06             	shl    $0x6,%eax
   1369e:	41 89 c0             	mov    %eax,%r8d
   136a1:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   136a8:	48 8b 08             	mov    (%rax),%rcx
            wholeNParticles*sizeof(particle_t)/gridX,
   136ab:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   136b1:	48 98                	cltq   
   136b3:	48 c1 e0 06          	shl    $0x6,%rax
   136b7:	48 89 c7             	mov    %rax,%rdi
   136ba:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   136be:	8b 00                	mov    (%rax),%eax
   136c0:	4c 63 d8             	movslq %eax,%r11
   136c3:	48 89 f8             	mov    %rdi,%rax
   136c6:	ba 00 00 00 00       	mov    $0x0,%edx
   136cb:	49 f7 f3             	div    %r11
        MPI_Scatter(
   136ce:	89 c7                	mov    %eax,%edi
   136d0:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
   136d7:	56                   	push   %rsi
   136d8:	6a 00                	push   $0x0
   136da:	4c 8d 0d 3f a1 00 00 	lea    0xa13f(%rip),%r9        # 1d820 <ompi_mpi_byte@@Base>
   136e1:	48 8d 15 38 a1 00 00 	lea    0xa138(%rip),%rdx        # 1d820 <ompi_mpi_byte@@Base>
   136e8:	89 fe                	mov    %edi,%esi
   136ea:	48 89 c7             	mov    %rax,%rdi
   136ed:	e8 7e 6a ff ff       	call   a170 <MPI_Scatter@plt>
   136f2:	48 83 c4 10          	add    $0x10,%rsp
            MPI_BYTE,
            0,
            topology->scatterXComm
        );

        #pragma omp parallel for
   136f6:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   136fa:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   136fe:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   13702:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   13706:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   1370a:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   1370e:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   13712:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   13716:	8b 85 74 ff ff ff    	mov    -0x8c(%rbp),%eax
   1371c:	89 45 f4             	mov    %eax,-0xc(%rbp)
   1371f:	8b 85 70 ff ff ff    	mov    -0x90(%rbp),%eax
   13725:	89 45 f0             	mov    %eax,-0x10(%rbp)
   13728:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   1372f:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   13733:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   1373a:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   1373e:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   13742:	b9 00 00 00 00       	mov    $0x0,%ecx
   13747:	ba 00 00 00 00       	mov    $0x0,%edx
   1374c:	48 89 c6             	mov    %rax,%rsi
   1374f:	48 8d 05 67 1a 00 00 	lea    0x1a67(%rip),%rax        # 151bd <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.0>
   13756:	48 89 c7             	mov    %rax,%rdi
   13759:	e8 52 70 ff ff       	call   a7b0 <GOMP_parallel@plt>
            ((*localParticlesX)+n)->velocity = vel+n*ndim;
            ((*localParticlesX)+n)->acceleration = acc+n*ndim;
            ((*localParticlesX)+n)->features = feat+n*nfeat;
        }

        MPI_Barrier(topology->scatterXComm);
   1375e:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13765:	48 8b 40 20          	mov    0x20(%rax),%rax
   13769:	48 89 c7             	mov    %rax,%rdi
   1376c:	e8 6f 6a ff ff       	call   a1e0 <MPI_Barrier@plt>
        MPI_Scatter(
   13771:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13778:	48 8b 48 20          	mov    0x20(%rax),%rcx
            wholeParticles->position,
            wholeNParticles*ndim/gridX,
            MPI_DOUBLE,
            pos, 
            (*localNParticlesX)*ndim,
   1377c:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   13783:	8b 00                	mov    (%rax),%eax
        MPI_Scatter(
   13785:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   1378c:	89 c6                	mov    %eax,%esi
            wholeNParticles*ndim/gridX,
   1378e:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   13794:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   1379b:	48 8b 55 88          	mov    -0x78(%rbp),%rdx
   1379f:	8b 3a                	mov    (%rdx),%edi
        MPI_Scatter(
   137a1:	99                   	cltd   
   137a2:	f7 ff                	idiv   %edi
   137a4:	89 c7                	mov    %eax,%edi
   137a6:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
            wholeParticles->position,
   137ad:	48 8b 40 08          	mov    0x8(%rax),%rax
        MPI_Scatter(
   137b1:	48 8b 55 a0          	mov    -0x60(%rbp),%rdx
   137b5:	51                   	push   %rcx
   137b6:	6a 00                	push   $0x0
   137b8:	4c 8d 0d 61 a4 00 00 	lea    0xa461(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   137bf:	41 89 f0             	mov    %esi,%r8d
   137c2:	48 89 d1             	mov    %rdx,%rcx
   137c5:	48 8d 15 54 a4 00 00 	lea    0xa454(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   137cc:	89 fe                	mov    %edi,%esi
   137ce:	48 89 c7             	mov    %rax,%rdi
   137d1:	e8 9a 69 ff ff       	call   a170 <MPI_Scatter@plt>
   137d6:	48 83 c4 10          	add    $0x10,%rsp
            MPI_DOUBLE,
            0,
            topology->scatterXComm
        );

        MPI_Barrier(topology->scatterXComm);
   137da:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   137e1:	48 8b 40 20          	mov    0x20(%rax),%rax
   137e5:	48 89 c7             	mov    %rax,%rdi
   137e8:	e8 f3 69 ff ff       	call   a1e0 <MPI_Barrier@plt>
        printf("(%d, %d) %s %d %d@%d\n", topology->rankX, topology->rankY, __FILE__, __LINE__, *localNParticlesX, wholeNParticles);
   137ed:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   137f4:	8b 30                	mov    (%rax),%esi
   137f6:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   137fd:	8b 50 04             	mov    0x4(%rax),%edx
   13800:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13807:	8b 00                	mov    (%rax),%eax
   13809:	48 83 ec 08          	sub    $0x8,%rsp
   1380d:	8b 8d 34 ff ff ff    	mov    -0xcc(%rbp),%ecx
   13813:	51                   	push   %rcx
   13814:	41 89 f1             	mov    %esi,%r9d
   13817:	41 b8 17 01 00 00    	mov    $0x117,%r8d
   1381d:	48 8d 0d b4 34 00 00 	lea    0x34b4(%rip),%rcx        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   13824:	89 c6                	mov    %eax,%esi
   13826:	48 8d 05 95 36 00 00 	lea    0x3695(%rip),%rax        # 16ec2 <_ZN3MPIL11LOCK_SHAREDE+0x222>
   1382d:	48 89 c7             	mov    %rax,%rdi
   13830:	b8 00 00 00 00       	mov    $0x0,%eax
   13835:	e8 a6 65 ff ff       	call   9de0 <printf@plt>
   1383a:	48 83 c4 10          	add    $0x10,%rsp
        MPI_Scatter(
   1383e:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13845:	48 8b 70 20          	mov    0x20(%rax),%rsi
            wholeParticles->velocity,
            wholeNParticles*ndim/gridX,
            MPI_DOUBLE,
            (*localParticlesX)->velocity, 
            (*localNParticlesX)*ndim,
   13849:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   13850:	8b 00                	mov    (%rax),%eax
        MPI_Scatter(
   13852:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   13859:	89 c7                	mov    %eax,%edi
            (*localParticlesX)->velocity, 
   1385b:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   13862:	48 8b 00             	mov    (%rax),%rax
   13865:	48 8b 48 10          	mov    0x10(%rax),%rcx
            wholeNParticles*ndim/gridX,
   13869:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   1386f:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   13876:	48 8b 55 88          	mov    -0x78(%rbp),%rdx
   1387a:	44 8b 1a             	mov    (%rdx),%r11d
        MPI_Scatter(
   1387d:	99                   	cltd   
   1387e:	41 f7 fb             	idiv   %r11d
   13881:	41 89 c2             	mov    %eax,%r10d
   13884:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
            wholeParticles->velocity,
   1388b:	48 8b 40 10          	mov    0x10(%rax),%rax
        MPI_Scatter(
   1388f:	56                   	push   %rsi
   13890:	6a 00                	push   $0x0
   13892:	4c 8d 0d 87 a3 00 00 	lea    0xa387(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   13899:	41 89 f8             	mov    %edi,%r8d
   1389c:	48 8d 15 7d a3 00 00 	lea    0xa37d(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   138a3:	44 89 d6             	mov    %r10d,%esi
   138a6:	48 89 c7             	mov    %rax,%rdi
   138a9:	e8 c2 68 ff ff       	call   a170 <MPI_Scatter@plt>
   138ae:	48 83 c4 10          	add    $0x10,%rsp
            MPI_DOUBLE,
            0,
            topology->scatterXComm
        );
        
        MPI_Barrier(topology->scatterXComm);
   138b2:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   138b9:	48 8b 40 20          	mov    0x20(%rax),%rax
   138bd:	48 89 c7             	mov    %rax,%rdi
   138c0:	e8 1b 69 ff ff       	call   a1e0 <MPI_Barrier@plt>
        printf("(%d, %d) %s %d %d@%d\n", topology->rankX, topology->rankY, __FILE__, __LINE__, *localNParticlesX, wholeNParticles);
   138c5:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   138cc:	8b 30                	mov    (%rax),%esi
   138ce:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   138d5:	8b 50 04             	mov    0x4(%rax),%edx
   138d8:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   138df:	8b 00                	mov    (%rax),%eax
   138e1:	48 83 ec 08          	sub    $0x8,%rsp
   138e5:	8b 8d 34 ff ff ff    	mov    -0xcc(%rbp),%ecx
   138eb:	51                   	push   %rcx
   138ec:	41 89 f1             	mov    %esi,%r9d
   138ef:	41 b8 24 01 00 00    	mov    $0x124,%r8d
   138f5:	48 8d 0d dc 33 00 00 	lea    0x33dc(%rip),%rcx        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   138fc:	89 c6                	mov    %eax,%esi
   138fe:	48 8d 05 bd 35 00 00 	lea    0x35bd(%rip),%rax        # 16ec2 <_ZN3MPIL11LOCK_SHAREDE+0x222>
   13905:	48 89 c7             	mov    %rax,%rdi
   13908:	b8 00 00 00 00       	mov    $0x0,%eax
   1390d:	e8 ce 64 ff ff       	call   9de0 <printf@plt>
   13912:	48 83 c4 10          	add    $0x10,%rsp
        MPI_Scatter(
   13916:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   1391d:	48 8b 70 20          	mov    0x20(%rax),%rsi
            wholeParticles->acceleration,
            wholeNParticles*ndim/gridX,
            MPI_DOUBLE,
            (*localParticlesX)->acceleration, 
            (*localNParticlesX)*ndim,
   13921:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   13928:	8b 00                	mov    (%rax),%eax
        MPI_Scatter(
   1392a:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   13931:	89 c7                	mov    %eax,%edi
            (*localParticlesX)->acceleration, 
   13933:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   1393a:	48 8b 00             	mov    (%rax),%rax
   1393d:	48 8b 48 18          	mov    0x18(%rax),%rcx
            wholeNParticles*ndim/gridX,
   13941:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   13947:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   1394e:	48 8b 55 88          	mov    -0x78(%rbp),%rdx
   13952:	44 8b 1a             	mov    (%rdx),%r11d
        MPI_Scatter(
   13955:	99                   	cltd   
   13956:	41 f7 fb             	idiv   %r11d
   13959:	41 89 c2             	mov    %eax,%r10d
   1395c:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
            wholeParticles->acceleration,
   13963:	48 8b 40 18          	mov    0x18(%rax),%rax
        MPI_Scatter(
   13967:	56                   	push   %rsi
   13968:	6a 00                	push   $0x0
   1396a:	4c 8d 0d af a2 00 00 	lea    0xa2af(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   13971:	41 89 f8             	mov    %edi,%r8d
   13974:	48 8d 15 a5 a2 00 00 	lea    0xa2a5(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   1397b:	44 89 d6             	mov    %r10d,%esi
   1397e:	48 89 c7             	mov    %rax,%rdi
   13981:	e8 ea 67 ff ff       	call   a170 <MPI_Scatter@plt>
   13986:	48 83 c4 10          	add    $0x10,%rsp
            MPI_DOUBLE,
            0,
            topology->scatterXComm
        );
        
        MPI_Barrier(topology->scatterXComm);
   1398a:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13991:	48 8b 40 20          	mov    0x20(%rax),%rax
   13995:	48 89 c7             	mov    %rax,%rdi
   13998:	e8 43 68 ff ff       	call   a1e0 <MPI_Barrier@plt>
        printf("(%d, %d) %s %d %d@%d\n", topology->rankX, topology->rankY, __FILE__, __LINE__, *localNParticlesX, wholeNParticles);
   1399d:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   139a4:	8b 30                	mov    (%rax),%esi
   139a6:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   139ad:	8b 50 04             	mov    0x4(%rax),%edx
   139b0:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   139b7:	8b 00                	mov    (%rax),%eax
   139b9:	48 83 ec 08          	sub    $0x8,%rsp
   139bd:	8b 8d 34 ff ff ff    	mov    -0xcc(%rbp),%ecx
   139c3:	51                   	push   %rcx
   139c4:	41 89 f1             	mov    %esi,%r9d
   139c7:	41 b8 31 01 00 00    	mov    $0x131,%r8d
   139cd:	48 8d 0d 04 33 00 00 	lea    0x3304(%rip),%rcx        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   139d4:	89 c6                	mov    %eax,%esi
   139d6:	48 8d 05 e5 34 00 00 	lea    0x34e5(%rip),%rax        # 16ec2 <_ZN3MPIL11LOCK_SHAREDE+0x222>
   139dd:	48 89 c7             	mov    %rax,%rdi
   139e0:	b8 00 00 00 00       	mov    $0x0,%eax
   139e5:	e8 f6 63 ff ff       	call   9de0 <printf@plt>
   139ea:	48 83 c4 10          	add    $0x10,%rsp
        MPI_Scatter(
   139ee:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   139f5:	48 8b 70 20          	mov    0x20(%rax),%rsi
            wholeParticles->features,
            wholeNParticles*nfeat/gridX,
            MPI_DOUBLE,
            (*localParticlesX)->features, 
            (*localNParticlesX)*nfeat,
   139f9:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   13a00:	8b 00                	mov    (%rax),%eax
        MPI_Scatter(
   13a02:	0f af 85 74 ff ff ff 	imul   -0x8c(%rbp),%eax
   13a09:	89 c7                	mov    %eax,%edi
            (*localParticlesX)->features, 
   13a0b:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   13a12:	48 8b 00             	mov    (%rax),%rax
   13a15:	48 8b 48 20          	mov    0x20(%rax),%rcx
            wholeNParticles*nfeat/gridX,
   13a19:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   13a1f:	0f af 85 74 ff ff ff 	imul   -0x8c(%rbp),%eax
   13a26:	48 8b 55 88          	mov    -0x78(%rbp),%rdx
   13a2a:	44 8b 1a             	mov    (%rdx),%r11d
        MPI_Scatter(
   13a2d:	99                   	cltd   
   13a2e:	41 f7 fb             	idiv   %r11d
   13a31:	41 89 c2             	mov    %eax,%r10d
   13a34:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
            wholeParticles->features,
   13a3b:	48 8b 40 20          	mov    0x20(%rax),%rax
        MPI_Scatter(
   13a3f:	56                   	push   %rsi
   13a40:	6a 00                	push   $0x0
   13a42:	4c 8d 0d d7 a1 00 00 	lea    0xa1d7(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   13a49:	41 89 f8             	mov    %edi,%r8d
   13a4c:	48 8d 15 cd a1 00 00 	lea    0xa1cd(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13a53:	44 89 d6             	mov    %r10d,%esi
   13a56:	48 89 c7             	mov    %rax,%rdi
   13a59:	e8 12 67 ff ff       	call   a170 <MPI_Scatter@plt>
   13a5e:	48 83 c4 10          	add    $0x10,%rsp
            0,
            topology->scatterXComm
        );
    }

    pos = (*localParticlesX)->position;
   13a62:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   13a69:	48 8b 00             	mov    (%rax),%rax
   13a6c:	48 8b 40 08          	mov    0x8(%rax),%rax
   13a70:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
    vel = (*localParticlesX)->velocity;
   13a74:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   13a7b:	48 8b 00             	mov    (%rax),%rax
   13a7e:	48 8b 40 10          	mov    0x10(%rax),%rax
   13a82:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
    acc = (*localParticlesX)->acceleration;
   13a86:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   13a8d:	48 8b 00             	mov    (%rax),%rax
   13a90:	48 8b 40 18          	mov    0x18(%rax),%rax
   13a94:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
    feat = (*localParticlesX)->features;
   13a98:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   13a9f:	48 8b 00             	mov    (%rax),%rax
   13aa2:	48 8b 40 20          	mov    0x20(%rax),%rax
   13aa6:	48 89 45 b8          	mov    %rax,-0x48(%rbp)

    MPI_Bcast(
   13aaa:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13ab1:	48 8b 50 10          	mov    0x10(%rax),%rdx
        (*localParticlesX),
        (*localNParticlesX)*sizeof(particle_t),
   13ab5:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   13abc:	8b 00                	mov    (%rax),%eax
   13abe:	48 98                	cltq   
    MPI_Bcast(
   13ac0:	c1 e0 06             	shl    $0x6,%eax
   13ac3:	89 c6                	mov    %eax,%esi
   13ac5:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   13acc:	48 8b 00             	mov    (%rax),%rax
   13acf:	49 89 d0             	mov    %rdx,%r8
   13ad2:	b9 00 00 00 00       	mov    $0x0,%ecx
   13ad7:	48 8d 15 42 9d 00 00 	lea    0x9d42(%rip),%rdx        # 1d820 <ompi_mpi_byte@@Base>
   13ade:	48 89 c7             	mov    %rax,%rdi
   13ae1:	e8 8a 63 ff ff       	call   9e70 <MPI_Bcast@plt>
        MPI_BYTE,
        0,
        topology->reduceXComm
    );

    #pragma omp parallel for
   13ae6:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   13aea:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   13aee:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   13af2:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   13af6:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   13afa:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   13afe:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   13b02:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   13b06:	8b 85 74 ff ff ff    	mov    -0x8c(%rbp),%eax
   13b0c:	89 45 f4             	mov    %eax,-0xc(%rbp)
   13b0f:	8b 85 70 ff ff ff    	mov    -0x90(%rbp),%eax
   13b15:	89 45 f0             	mov    %eax,-0x10(%rbp)
   13b18:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   13b1f:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   13b23:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   13b2a:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   13b2e:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   13b32:	b9 00 00 00 00       	mov    $0x0,%ecx
   13b37:	ba 00 00 00 00       	mov    $0x0,%edx
   13b3c:	48 89 c6             	mov    %rax,%rsi
   13b3f:	48 8d 05 09 18 00 00 	lea    0x1809(%rip),%rax        # 1534f <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.1>
   13b46:	48 89 c7             	mov    %rax,%rdi
   13b49:	e8 62 6c ff ff       	call   a7b0 <GOMP_parallel@plt>
        (*localParticlesX+n)->position = pos+n*ndim;
        (*localParticlesX+n)->velocity = vel+n*ndim;
        (*localParticlesX+n)->acceleration = acc+n*ndim;
        (*localParticlesX+n)->features = feat+n*nfeat;
    }
    MPI_Bcast(
   13b4e:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13b55:	48 8b 50 10          	mov    0x10(%rax),%rdx
        (*localParticlesX)->position,
        (*localNParticlesX)*ndim,
   13b59:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   13b60:	8b 00                	mov    (%rax),%eax
    MPI_Bcast(
   13b62:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   13b69:	89 c6                	mov    %eax,%esi
        (*localParticlesX)->position,
   13b6b:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   13b72:	48 8b 00             	mov    (%rax),%rax
   13b75:	48 8b 40 08          	mov    0x8(%rax),%rax
    MPI_Bcast(
   13b79:	49 89 d0             	mov    %rdx,%r8
   13b7c:	b9 00 00 00 00       	mov    $0x0,%ecx
   13b81:	48 8d 15 98 a0 00 00 	lea    0xa098(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13b88:	48 89 c7             	mov    %rax,%rdi
   13b8b:	e8 e0 62 ff ff       	call   9e70 <MPI_Bcast@plt>
        MPI_DOUBLE,
        0,
        topology->reduceXComm
    );
    MPI_Bcast(
   13b90:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13b97:	48 8b 50 10          	mov    0x10(%rax),%rdx
        (*localParticlesX)->velocity,
        (*localNParticlesX)*ndim,
   13b9b:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   13ba2:	8b 00                	mov    (%rax),%eax
    MPI_Bcast(
   13ba4:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   13bab:	89 c6                	mov    %eax,%esi
        (*localParticlesX)->velocity,
   13bad:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   13bb4:	48 8b 00             	mov    (%rax),%rax
   13bb7:	48 8b 40 10          	mov    0x10(%rax),%rax
    MPI_Bcast(
   13bbb:	49 89 d0             	mov    %rdx,%r8
   13bbe:	b9 00 00 00 00       	mov    $0x0,%ecx
   13bc3:	48 8d 15 56 a0 00 00 	lea    0xa056(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13bca:	48 89 c7             	mov    %rax,%rdi
   13bcd:	e8 9e 62 ff ff       	call   9e70 <MPI_Bcast@plt>
        MPI_DOUBLE,
        0,
        topology->reduceXComm
    );
    MPI_Bcast(
   13bd2:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13bd9:	48 8b 50 10          	mov    0x10(%rax),%rdx
        (*localParticlesX)->acceleration,
        (*localNParticlesX)*ndim,
   13bdd:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   13be4:	8b 00                	mov    (%rax),%eax
    MPI_Bcast(
   13be6:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   13bed:	89 c6                	mov    %eax,%esi
        (*localParticlesX)->acceleration,
   13bef:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   13bf6:	48 8b 00             	mov    (%rax),%rax
   13bf9:	48 8b 40 18          	mov    0x18(%rax),%rax
    MPI_Bcast(
   13bfd:	49 89 d0             	mov    %rdx,%r8
   13c00:	b9 00 00 00 00       	mov    $0x0,%ecx
   13c05:	48 8d 15 14 a0 00 00 	lea    0xa014(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13c0c:	48 89 c7             	mov    %rax,%rdi
   13c0f:	e8 5c 62 ff ff       	call   9e70 <MPI_Bcast@plt>
        MPI_DOUBLE,
        0,
        topology->reduceXComm
    );
    MPI_Bcast(
   13c14:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13c1b:	48 8b 50 10          	mov    0x10(%rax),%rdx
        (*localParticlesX)->features,
        (*localNParticlesX)*nfeat,
   13c1f:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   13c26:	8b 00                	mov    (%rax),%eax
    MPI_Bcast(
   13c28:	0f af 85 74 ff ff ff 	imul   -0x8c(%rbp),%eax
   13c2f:	89 c6                	mov    %eax,%esi
        (*localParticlesX)->features,
   13c31:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   13c38:	48 8b 00             	mov    (%rax),%rax
   13c3b:	48 8b 40 20          	mov    0x20(%rax),%rax
    MPI_Bcast(
   13c3f:	49 89 d0             	mov    %rdx,%r8
   13c42:	b9 00 00 00 00       	mov    $0x0,%ecx
   13c47:	48 8d 15 d2 9f 00 00 	lea    0x9fd2(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13c4e:	48 89 c7             	mov    %rax,%rdi
   13c51:	e8 1a 62 ff ff       	call   9e70 <MPI_Bcast@plt>
        0,
        topology->reduceXComm
    );


    if(rankX==0) {
   13c56:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13c5d:	8b 00                	mov    (%rax),%eax
   13c5f:	85 c0                	test   %eax,%eax
   13c61:	0f 85 f1 02 00 00    	jne    13f58 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0xc8b>
        pos = (*localParticlesY)->position;
   13c67:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13c6e:	48 8b 00             	mov    (%rax),%rax
   13c71:	48 8b 40 08          	mov    0x8(%rax),%rax
   13c75:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
        vel = (*localParticlesY)->velocity;
   13c79:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13c80:	48 8b 00             	mov    (%rax),%rax
   13c83:	48 8b 40 10          	mov    0x10(%rax),%rax
   13c87:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
        acc = (*localParticlesY)->acceleration;
   13c8b:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13c92:	48 8b 00             	mov    (%rax),%rax
   13c95:	48 8b 40 18          	mov    0x18(%rax),%rax
   13c99:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
        feat = (*localParticlesY)->features;
   13c9d:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13ca4:	48 8b 00             	mov    (%rax),%rax
   13ca7:	48 8b 40 20          	mov    0x20(%rax),%rax
   13cab:	48 89 45 b8          	mov    %rax,-0x48(%rbp)

        MPI_Scatter(
   13caf:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13cb6:	48 8b 70 28          	mov    0x28(%rax),%rsi
            wholeParticles,
            wholeNParticles*sizeof(particle_t)/gridY,
            MPI_BYTE,
            (*localParticlesY),
            (*localNParticlesY)*sizeof(particle_t),
   13cba:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   13cc1:	8b 00                	mov    (%rax),%eax
   13cc3:	48 98                	cltq   
        MPI_Scatter(
   13cc5:	c1 e0 06             	shl    $0x6,%eax
   13cc8:	41 89 c0             	mov    %eax,%r8d
   13ccb:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13cd2:	48 8b 08             	mov    (%rax),%rcx
            wholeNParticles*sizeof(particle_t)/gridY,
   13cd5:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   13cdb:	48 98                	cltq   
   13cdd:	48 c1 e0 06          	shl    $0x6,%rax
   13ce1:	48 89 c7             	mov    %rax,%rdi
   13ce4:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   13ce8:	8b 00                	mov    (%rax),%eax
   13cea:	4c 63 d8             	movslq %eax,%r11
   13ced:	48 89 f8             	mov    %rdi,%rax
   13cf0:	ba 00 00 00 00       	mov    $0x0,%edx
   13cf5:	49 f7 f3             	div    %r11
        MPI_Scatter(
   13cf8:	89 c7                	mov    %eax,%edi
   13cfa:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
   13d01:	56                   	push   %rsi
   13d02:	6a 00                	push   $0x0
   13d04:	4c 8d 0d 15 9b 00 00 	lea    0x9b15(%rip),%r9        # 1d820 <ompi_mpi_byte@@Base>
   13d0b:	48 8d 15 0e 9b 00 00 	lea    0x9b0e(%rip),%rdx        # 1d820 <ompi_mpi_byte@@Base>
   13d12:	89 fe                	mov    %edi,%esi
   13d14:	48 89 c7             	mov    %rax,%rdi
   13d17:	e8 54 64 ff ff       	call   a170 <MPI_Scatter@plt>
   13d1c:	48 83 c4 10          	add    $0x10,%rsp
            MPI_BYTE,
            0,
            topology->scatterYComm
        );

        #pragma omp parallel for
   13d20:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   13d24:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   13d28:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   13d2c:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   13d30:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   13d34:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   13d38:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   13d3c:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   13d40:	8b 85 74 ff ff ff    	mov    -0x8c(%rbp),%eax
   13d46:	89 45 f4             	mov    %eax,-0xc(%rbp)
   13d49:	8b 85 70 ff ff ff    	mov    -0x90(%rbp),%eax
   13d4f:	89 45 f0             	mov    %eax,-0x10(%rbp)
   13d52:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   13d59:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   13d5d:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13d64:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   13d68:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   13d6c:	b9 00 00 00 00       	mov    $0x0,%ecx
   13d71:	ba 00 00 00 00       	mov    $0x0,%edx
   13d76:	48 89 c6             	mov    %rax,%rsi
   13d79:	48 8d 05 61 17 00 00 	lea    0x1761(%rip),%rax        # 154e1 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.2>
   13d80:	48 89 c7             	mov    %rax,%rdi
   13d83:	e8 28 6a ff ff       	call   a7b0 <GOMP_parallel@plt>
            (*localParticlesY+n)->velocity = vel+n*ndim;
            (*localParticlesY+n)->acceleration = acc+n*ndim;
            (*localParticlesY+n)->features = feat+n*nfeat;
        }

        MPI_Scatter(
   13d88:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13d8f:	48 8b 70 28          	mov    0x28(%rax),%rsi
            wholeParticles->position,
            wholeNParticles*ndim/gridY,
            MPI_DOUBLE,
            (*localParticlesY)->position, 
            (*localNParticlesY)*ndim,
   13d93:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   13d9a:	8b 00                	mov    (%rax),%eax
        MPI_Scatter(
   13d9c:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   13da3:	89 c7                	mov    %eax,%edi
            (*localParticlesY)->position, 
   13da5:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13dac:	48 8b 00             	mov    (%rax),%rax
   13daf:	48 8b 48 08          	mov    0x8(%rax),%rcx
            wholeNParticles*ndim/gridY,
   13db3:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   13db9:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   13dc0:	48 8b 55 90          	mov    -0x70(%rbp),%rdx
   13dc4:	44 8b 1a             	mov    (%rdx),%r11d
        MPI_Scatter(
   13dc7:	99                   	cltd   
   13dc8:	41 f7 fb             	idiv   %r11d
   13dcb:	41 89 c2             	mov    %eax,%r10d
   13dce:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
            wholeParticles->position,
   13dd5:	48 8b 40 08          	mov    0x8(%rax),%rax
        MPI_Scatter(
   13dd9:	56                   	push   %rsi
   13dda:	6a 00                	push   $0x0
   13ddc:	4c 8d 0d 3d 9e 00 00 	lea    0x9e3d(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   13de3:	41 89 f8             	mov    %edi,%r8d
   13de6:	48 8d 15 33 9e 00 00 	lea    0x9e33(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13ded:	44 89 d6             	mov    %r10d,%esi
   13df0:	48 89 c7             	mov    %rax,%rdi
   13df3:	e8 78 63 ff ff       	call   a170 <MPI_Scatter@plt>
   13df8:	48 83 c4 10          	add    $0x10,%rsp
            MPI_DOUBLE,
            0,
            topology->scatterYComm
        );
        MPI_Scatter(
   13dfc:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13e03:	48 8b 70 28          	mov    0x28(%rax),%rsi
            wholeParticles->velocity,
            wholeNParticles*ndim/gridY,
            MPI_DOUBLE,
            (*localParticlesY)->velocity, 
            (*localNParticlesY)*ndim,
   13e07:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   13e0e:	8b 00                	mov    (%rax),%eax
        MPI_Scatter(
   13e10:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   13e17:	89 c7                	mov    %eax,%edi
            (*localParticlesY)->velocity, 
   13e19:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13e20:	48 8b 00             	mov    (%rax),%rax
   13e23:	48 8b 48 10          	mov    0x10(%rax),%rcx
            wholeNParticles*ndim/gridY,
   13e27:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   13e2d:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   13e34:	48 8b 55 90          	mov    -0x70(%rbp),%rdx
   13e38:	44 8b 1a             	mov    (%rdx),%r11d
        MPI_Scatter(
   13e3b:	99                   	cltd   
   13e3c:	41 f7 fb             	idiv   %r11d
   13e3f:	41 89 c2             	mov    %eax,%r10d
   13e42:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
            wholeParticles->velocity,
   13e49:	48 8b 40 10          	mov    0x10(%rax),%rax
        MPI_Scatter(
   13e4d:	56                   	push   %rsi
   13e4e:	6a 00                	push   $0x0
   13e50:	4c 8d 0d c9 9d 00 00 	lea    0x9dc9(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   13e57:	41 89 f8             	mov    %edi,%r8d
   13e5a:	48 8d 15 bf 9d 00 00 	lea    0x9dbf(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13e61:	44 89 d6             	mov    %r10d,%esi
   13e64:	48 89 c7             	mov    %rax,%rdi
   13e67:	e8 04 63 ff ff       	call   a170 <MPI_Scatter@plt>
   13e6c:	48 83 c4 10          	add    $0x10,%rsp
            MPI_DOUBLE,
            0,
            topology->scatterYComm
        );
        MPI_Scatter(
   13e70:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13e77:	48 8b 70 28          	mov    0x28(%rax),%rsi
            wholeParticles->acceleration,
            wholeNParticles*ndim/gridY,
            MPI_DOUBLE,
            (*localParticlesY)->acceleration, 
            (*localNParticlesY)*ndim,
   13e7b:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   13e82:	8b 00                	mov    (%rax),%eax
        MPI_Scatter(
   13e84:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   13e8b:	89 c7                	mov    %eax,%edi
            (*localParticlesY)->acceleration, 
   13e8d:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13e94:	48 8b 00             	mov    (%rax),%rax
   13e97:	48 8b 48 18          	mov    0x18(%rax),%rcx
            wholeNParticles*ndim/gridY,
   13e9b:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   13ea1:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   13ea8:	48 8b 55 90          	mov    -0x70(%rbp),%rdx
   13eac:	44 8b 12             	mov    (%rdx),%r10d
        MPI_Scatter(
   13eaf:	99                   	cltd   
   13eb0:	41 f7 fa             	idiv   %r10d
   13eb3:	41 89 c2             	mov    %eax,%r10d
   13eb6:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
            wholeParticles->acceleration,
   13ebd:	48 8b 40 18          	mov    0x18(%rax),%rax
        MPI_Scatter(
   13ec1:	56                   	push   %rsi
   13ec2:	6a 00                	push   $0x0
   13ec4:	4c 8d 0d 55 9d 00 00 	lea    0x9d55(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   13ecb:	41 89 f8             	mov    %edi,%r8d
   13ece:	48 8d 15 4b 9d 00 00 	lea    0x9d4b(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13ed5:	44 89 d6             	mov    %r10d,%esi
   13ed8:	48 89 c7             	mov    %rax,%rdi
   13edb:	e8 90 62 ff ff       	call   a170 <MPI_Scatter@plt>
   13ee0:	48 83 c4 10          	add    $0x10,%rsp
            MPI_DOUBLE,
            0,
            topology->scatterYComm
        );
        MPI_Scatter(
   13ee4:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13eeb:	48 8b 70 28          	mov    0x28(%rax),%rsi
            wholeParticles->features,
            wholeNParticles*nfeat/gridY,
            MPI_DOUBLE,
            (*localParticlesY)->features, 
            (*localNParticlesY)*nfeat,
   13eef:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   13ef6:	8b 00                	mov    (%rax),%eax
        MPI_Scatter(
   13ef8:	0f af 85 74 ff ff ff 	imul   -0x8c(%rbp),%eax
   13eff:	89 c7                	mov    %eax,%edi
            (*localParticlesY)->features, 
   13f01:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13f08:	48 8b 00             	mov    (%rax),%rax
   13f0b:	48 8b 48 20          	mov    0x20(%rax),%rcx
            wholeNParticles*nfeat/gridY,
   13f0f:	8b 85 34 ff ff ff    	mov    -0xcc(%rbp),%eax
   13f15:	0f af 85 74 ff ff ff 	imul   -0x8c(%rbp),%eax
   13f1c:	48 8b 55 90          	mov    -0x70(%rbp),%rdx
   13f20:	44 8b 1a             	mov    (%rdx),%r11d
        MPI_Scatter(
   13f23:	99                   	cltd   
   13f24:	41 f7 fb             	idiv   %r11d
   13f27:	41 89 c2             	mov    %eax,%r10d
   13f2a:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
            wholeParticles->features,
   13f31:	48 8b 40 20          	mov    0x20(%rax),%rax
        MPI_Scatter(
   13f35:	56                   	push   %rsi
   13f36:	6a 00                	push   $0x0
   13f38:	4c 8d 0d e1 9c 00 00 	lea    0x9ce1(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   13f3f:	41 89 f8             	mov    %edi,%r8d
   13f42:	48 8d 15 d7 9c 00 00 	lea    0x9cd7(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13f49:	44 89 d6             	mov    %r10d,%esi
   13f4c:	48 89 c7             	mov    %rax,%rdi
   13f4f:	e8 1c 62 ff ff       	call   a170 <MPI_Scatter@plt>
   13f54:	48 83 c4 10          	add    $0x10,%rsp
            0,
            topology->scatterYComm
        );
    }

    pos = (*localParticlesY)->position;
   13f58:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13f5f:	48 8b 00             	mov    (%rax),%rax
   13f62:	48 8b 40 08          	mov    0x8(%rax),%rax
   13f66:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
    vel = (*localParticlesY)->velocity;
   13f6a:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13f71:	48 8b 00             	mov    (%rax),%rax
   13f74:	48 8b 40 10          	mov    0x10(%rax),%rax
   13f78:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
    acc = (*localParticlesY)->acceleration;
   13f7c:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13f83:	48 8b 00             	mov    (%rax),%rax
   13f86:	48 8b 40 18          	mov    0x18(%rax),%rax
   13f8a:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
    feat = (*localParticlesY)->features;
   13f8e:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13f95:	48 8b 00             	mov    (%rax),%rax
   13f98:	48 8b 40 20          	mov    0x20(%rax),%rax
   13f9c:	48 89 45 b8          	mov    %rax,-0x48(%rbp)

    MPI_Bcast(
   13fa0:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   13fa7:	48 8b 50 18          	mov    0x18(%rax),%rdx
        (*localParticlesY),
        (*localNParticlesY)*sizeof(particle_t),
   13fab:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   13fb2:	8b 00                	mov    (%rax),%eax
   13fb4:	48 98                	cltq   
    MPI_Bcast(
   13fb6:	c1 e0 06             	shl    $0x6,%eax
   13fb9:	89 c6                	mov    %eax,%esi
   13fbb:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   13fc2:	48 8b 00             	mov    (%rax),%rax
   13fc5:	49 89 d0             	mov    %rdx,%r8
   13fc8:	b9 00 00 00 00       	mov    $0x0,%ecx
   13fcd:	48 8d 15 4c 98 00 00 	lea    0x984c(%rip),%rdx        # 1d820 <ompi_mpi_byte@@Base>
   13fd4:	48 89 c7             	mov    %rax,%rdi
   13fd7:	e8 94 5e ff ff       	call   9e70 <MPI_Bcast@plt>
        MPI_BYTE,
        0,
        topology->reduceYComm
    );

    #pragma omp parallel for
   13fdc:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   13fe0:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   13fe4:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   13fe8:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   13fec:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   13ff0:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   13ff4:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   13ff8:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   13ffc:	8b 85 74 ff ff ff    	mov    -0x8c(%rbp),%eax
   14002:	89 45 f4             	mov    %eax,-0xc(%rbp)
   14005:	8b 85 70 ff ff ff    	mov    -0x90(%rbp),%eax
   1400b:	89 45 f0             	mov    %eax,-0x10(%rbp)
   1400e:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   14015:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   14019:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   14020:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   14024:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   14028:	b9 00 00 00 00       	mov    $0x0,%ecx
   1402d:	ba 00 00 00 00       	mov    $0x0,%edx
   14032:	48 89 c6             	mov    %rax,%rsi
   14035:	48 8d 05 37 16 00 00 	lea    0x1637(%rip),%rax        # 15673 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.3>
   1403c:	48 89 c7             	mov    %rax,%rdi
   1403f:	e8 6c 67 ff ff       	call   a7b0 <GOMP_parallel@plt>
        (*localParticlesY+n)->position = pos+n*ndim;
        (*localParticlesY+n)->velocity = vel+n*ndim;
        (*localParticlesY+n)->acceleration = acc+n*ndim;
        (*localParticlesY+n)->features = feat+n*nfeat;
    }
    MPI_Bcast(
   14044:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   1404b:	48 8b 50 18          	mov    0x18(%rax),%rdx
        (*localParticlesY)->position,
        (*localNParticlesY)*ndim,
   1404f:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   14056:	8b 00                	mov    (%rax),%eax
    MPI_Bcast(
   14058:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   1405f:	89 c6                	mov    %eax,%esi
        (*localParticlesY)->position,
   14061:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   14068:	48 8b 00             	mov    (%rax),%rax
   1406b:	48 8b 40 08          	mov    0x8(%rax),%rax
    MPI_Bcast(
   1406f:	49 89 d0             	mov    %rdx,%r8
   14072:	b9 00 00 00 00       	mov    $0x0,%ecx
   14077:	48 8d 15 a2 9b 00 00 	lea    0x9ba2(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   1407e:	48 89 c7             	mov    %rax,%rdi
   14081:	e8 ea 5d ff ff       	call   9e70 <MPI_Bcast@plt>
        MPI_DOUBLE,
        0,
        topology->reduceYComm
    );
    MPI_Bcast(
   14086:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   1408d:	48 8b 50 18          	mov    0x18(%rax),%rdx
        (*localParticlesY)->velocity,
        (*localNParticlesY)*ndim,
   14091:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   14098:	8b 00                	mov    (%rax),%eax
    MPI_Bcast(
   1409a:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   140a1:	89 c6                	mov    %eax,%esi
        (*localParticlesY)->velocity,
   140a3:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   140aa:	48 8b 00             	mov    (%rax),%rax
   140ad:	48 8b 40 10          	mov    0x10(%rax),%rax
    MPI_Bcast(
   140b1:	49 89 d0             	mov    %rdx,%r8
   140b4:	b9 00 00 00 00       	mov    $0x0,%ecx
   140b9:	48 8d 15 60 9b 00 00 	lea    0x9b60(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   140c0:	48 89 c7             	mov    %rax,%rdi
   140c3:	e8 a8 5d ff ff       	call   9e70 <MPI_Bcast@plt>
        MPI_DOUBLE,
        0,
        topology->reduceYComm
    );
    MPI_Bcast(
   140c8:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   140cf:	48 8b 50 18          	mov    0x18(%rax),%rdx
        (*localParticlesY)->acceleration,
        (*localNParticlesY)*ndim,
   140d3:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   140da:	8b 00                	mov    (%rax),%eax
    MPI_Bcast(
   140dc:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   140e3:	89 c6                	mov    %eax,%esi
        (*localParticlesY)->acceleration,
   140e5:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   140ec:	48 8b 00             	mov    (%rax),%rax
   140ef:	48 8b 40 18          	mov    0x18(%rax),%rax
    MPI_Bcast(
   140f3:	49 89 d0             	mov    %rdx,%r8
   140f6:	b9 00 00 00 00       	mov    $0x0,%ecx
   140fb:	48 8d 15 1e 9b 00 00 	lea    0x9b1e(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   14102:	48 89 c7             	mov    %rax,%rdi
   14105:	e8 66 5d ff ff       	call   9e70 <MPI_Bcast@plt>
        MPI_DOUBLE,
        0,
        topology->reduceYComm
    );
    MPI_Bcast(
   1410a:	48 8b 85 28 ff ff ff 	mov    -0xd8(%rbp),%rax
   14111:	48 8b 50 18          	mov    0x18(%rax),%rdx
        (*localParticlesY)->features,
        (*localNParticlesY)*nfeat,
   14115:	48 8b 85 40 ff ff ff 	mov    -0xc0(%rbp),%rax
   1411c:	8b 00                	mov    (%rax),%eax
    MPI_Bcast(
   1411e:	0f af 85 74 ff ff ff 	imul   -0x8c(%rbp),%eax
   14125:	89 c6                	mov    %eax,%esi
        (*localParticlesY)->features,
   14127:	48 8b 85 48 ff ff ff 	mov    -0xb8(%rbp),%rax
   1412e:	48 8b 00             	mov    (%rax),%rax
   14131:	48 8b 40 20          	mov    0x20(%rax),%rax
    MPI_Bcast(
   14135:	49 89 d0             	mov    %rdx,%r8
   14138:	b9 00 00 00 00       	mov    $0x0,%ecx
   1413d:	48 8d 15 dc 9a 00 00 	lea    0x9adc(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   14144:	48 89 c7             	mov    %rax,%rdi
   14147:	e8 24 5d ff ff       	call   9e70 <MPI_Bcast@plt>
        MPI_DOUBLE,
        0,
        topology->reduceYComm
    );
    if(wholeParticlesInitNull) {
   1414c:	80 bd 6f ff ff ff 00 	cmpb   $0x0,-0x91(%rbp)
   14153:	74 1a                	je     1416f <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0xea2>
        free_particles(wholeParticles);
   14155:	48 8b 85 38 ff ff ff 	mov    -0xc8(%rbp),%rax
   1415c:	48 89 c7             	mov    %rax,%rdi
   1415f:	e8 c4 da ff ff       	call   11c28 <_Z14free_particlesP10particle_t>
        wholeParticles = nullptr;
   14164:	48 c7 85 38 ff ff ff 	movq   $0x0,-0xc8(%rbp)
   1416b:	00 00 00 00 
    }
}
   1416f:	90                   	nop
   14170:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14174:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   1417b:	00 00 
   1417d:	74 05                	je     14184 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0xeb7>
   1417f:	e8 5c 63 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   14184:	c9                   	leave  
   14185:	c3                   	ret    

0000000000014186 <_Z7scatterPP17chunk_particles_tS1_S0_P10topology_t>:

void scatter(chunk_particles_t** localChunkParticlesX, 
                chunk_particles_t** localChunkParticlesY,
                chunk_particles_t* wholeChunkParticles,
                topology_t* topology) {
   14186:	f3 0f 1e fa          	endbr64 
   1418a:	55                   	push   %rbp
   1418b:	48 89 e5             	mov    %rsp,%rbp
   1418e:	48 83 ec 50          	sub    $0x50,%rsp
   14192:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
   14196:	48 89 75 c0          	mov    %rsi,-0x40(%rbp)
   1419a:	48 89 55 b8          	mov    %rdx,-0x48(%rbp)
   1419e:	48 89 4d b0          	mov    %rcx,-0x50(%rbp)
   141a2:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   141a9:	00 00 
   141ab:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   141af:	31 c0                	xor    %eax,%eax
    *localChunkParticlesX = static_cast<chunk_particles_t*>(malloc(sizeof(chunk_particles_t)));
   141b1:	bf 10 00 00 00       	mov    $0x10,%edi
   141b6:	e8 75 64 ff ff       	call   a630 <malloc@plt>
   141bb:	48 89 c2             	mov    %rax,%rdx
   141be:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   141c2:	48 89 10             	mov    %rdx,(%rax)
    *localChunkParticlesY = static_cast<chunk_particles_t*>(malloc(sizeof(chunk_particles_t)));
   141c5:	bf 10 00 00 00       	mov    $0x10,%edi
   141ca:	e8 61 64 ff ff       	call   a630 <malloc@plt>
   141cf:	48 89 c2             	mov    %rax,%rdx
   141d2:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   141d6:	48 89 10             	mov    %rdx,(%rax)
    particle_t *localParticlesX=nullptr, *localParticlesY=nullptr, *wholeParticles=nullptr;
   141d9:	48 c7 45 e0 00 00 00 	movq   $0x0,-0x20(%rbp)
   141e0:	00 
   141e1:	48 c7 45 e8 00 00 00 	movq   $0x0,-0x18(%rbp)
   141e8:	00 
   141e9:	48 c7 45 f0 00 00 00 	movq   $0x0,-0x10(%rbp)
   141f0:	00 
    int localNParticlesX=-1, localNParticlesY=-1, wholeNParticles;
   141f1:	c7 45 d4 ff ff ff ff 	movl   $0xffffffff,-0x2c(%rbp)
   141f8:	c7 45 d8 ff ff ff ff 	movl   $0xffffffff,-0x28(%rbp)
    if(wholeChunkParticles==nullptr){
   141ff:	48 83 7d b8 00       	cmpq   $0x0,-0x48(%rbp)
   14204:	75 11                	jne    14217 <_Z7scatterPP17chunk_particles_tS1_S0_P10topology_t+0x91>
        wholeParticles = nullptr;
   14206:	48 c7 45 f0 00 00 00 	movq   $0x0,-0x10(%rbp)
   1420d:	00 
        wholeNParticles = -1;
   1420e:	c7 45 dc ff ff ff ff 	movl   $0xffffffff,-0x24(%rbp)
   14215:	eb 15                	jmp    1422c <_Z7scatterPP17chunk_particles_tS1_S0_P10topology_t+0xa6>
    }
    else {
        wholeParticles = wholeChunkParticles->particles;
   14217:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   1421b:	48 8b 00             	mov    (%rax),%rax
   1421e:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
        wholeNParticles = wholeChunkParticles->nParticle;
   14222:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14226:	8b 40 08             	mov    0x8(%rax),%eax
   14229:	89 45 dc             	mov    %eax,-0x24(%rbp)
    }
    scatter(
   1422c:	44 8b 45 dc          	mov    -0x24(%rbp),%r8d
   14230:	48 8b 7d f0          	mov    -0x10(%rbp),%rdi
   14234:	48 8d 4d d8          	lea    -0x28(%rbp),%rcx
   14238:	48 8d 55 e8          	lea    -0x18(%rbp),%rdx
   1423c:	48 8d 75 d4          	lea    -0x2c(%rbp),%rsi
   14240:	48 8d 45 e0          	lea    -0x20(%rbp),%rax
   14244:	48 83 ec 08          	sub    $0x8,%rsp
   14248:	ff 75 b0             	push   -0x50(%rbp)
   1424b:	45 89 c1             	mov    %r8d,%r9d
   1424e:	49 89 f8             	mov    %rdi,%r8
   14251:	48 89 c7             	mov    %rax,%rdi
   14254:	e8 74 f0 ff ff       	call   132cd <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t>
   14259:	48 83 c4 10          	add    $0x10,%rsp
        &localParticlesX, &localNParticlesX, &localParticlesY, &localNParticlesY,
        wholeParticles, wholeNParticles, topology
    );
    (*localChunkParticlesX)->particles = localParticlesX;
   1425d:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14261:	48 8b 00             	mov    (%rax),%rax
   14264:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
   14268:	48 89 10             	mov    %rdx,(%rax)
    (*localChunkParticlesX)->nParticle = localNParticlesX;
   1426b:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   1426f:	48 8b 00             	mov    (%rax),%rax
   14272:	8b 55 d4             	mov    -0x2c(%rbp),%edx
   14275:	89 50 08             	mov    %edx,0x8(%rax)
    (*localChunkParticlesY)->particles = localParticlesY;
   14278:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   1427c:	48 8b 00             	mov    (%rax),%rax
   1427f:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   14283:	48 89 10             	mov    %rdx,(%rax)
    (*localChunkParticlesY)->nParticle = localNParticlesY;
   14286:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   1428a:	48 8b 00             	mov    (%rax),%rax
   1428d:	8b 55 d8             	mov    -0x28(%rbp),%edx
   14290:	89 50 08             	mov    %edx,0x8(%rax)
}
   14293:	90                   	nop
   14294:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14298:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   1429f:	00 00 
   142a1:	74 05                	je     142a8 <_Z7scatterPP17chunk_particles_tS1_S0_P10topology_t+0x122>
   142a3:	e8 38 62 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   142a8:	c9                   	leave  
   142a9:	c3                   	ret    

00000000000142aa <_Z6reduceP10particle_tiS0_iP10topology_t>:

void reduce(particle_t* localParticlesX, int localNParticlesX, 
        particle_t* localParticlesY, int localNParticlesY, topology_t* topology) {
   142aa:	f3 0f 1e fa          	endbr64 
   142ae:	55                   	push   %rbp
   142af:	48 89 e5             	mov    %rsp,%rbp
   142b2:	48 83 ec 50          	sub    $0x50,%rsp
   142b6:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
   142ba:	89 75 c4             	mov    %esi,-0x3c(%rbp)
   142bd:	48 89 55 b8          	mov    %rdx,-0x48(%rbp)
   142c1:	89 4d c0             	mov    %ecx,-0x40(%rbp)
   142c4:	4c 89 45 b0          	mov    %r8,-0x50(%rbp)
   142c8:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   142cf:	00 00 
   142d1:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   142d5:	31 c0                	xor    %eax,%eax
    double* recvBuf1;
    recvBuf1 = static_cast<double*>(malloc(localNParticlesX*localParticlesX->ndim*sizeof(double)));
   142d7:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   142db:	8b 40 30             	mov    0x30(%rax),%eax
   142de:	0f af 45 c4          	imul   -0x3c(%rbp),%eax
   142e2:	48 98                	cltq   
   142e4:	48 c1 e0 03          	shl    $0x3,%rax
   142e8:	48 89 c7             	mov    %rax,%rdi
   142eb:	e8 40 63 ff ff       	call   a630 <malloc@plt>
   142f0:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
    MPI_Allreduce(
   142f4:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   142f8:	48 8b 48 18          	mov    0x18(%rax),%rcx
        localParticlesX->acceleration,
        recvBuf1,
        localNParticlesX*localParticlesX->ndim,
   142fc:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14300:	8b 40 30             	mov    0x30(%rax),%eax
    MPI_Allreduce(
   14303:	0f af 45 c4          	imul   -0x3c(%rbp),%eax
   14307:	89 c2                	mov    %eax,%edx
        localParticlesX->acceleration,
   14309:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   1430d:	48 8b 40 18          	mov    0x18(%rax),%rax
    MPI_Allreduce(
   14311:	48 8b 75 d0          	mov    -0x30(%rbp),%rsi
   14315:	49 89 c9             	mov    %rcx,%r9
   14318:	4c 8d 05 01 8d 00 00 	lea    0x8d01(%rip),%r8        # 1d020 <ompi_mpi_op_sum@@Base>
   1431f:	48 8d 0d fa 98 00 00 	lea    0x98fa(%rip),%rcx        # 1dc20 <ompi_mpi_double@@Base>
   14326:	48 89 c7             	mov    %rax,%rdi
   14329:	e8 32 62 ff ff       	call   a560 <MPI_Allreduce@plt>
        MPI_DOUBLE,
        MPI_SUM,
        topology->reduceYComm
    );
    #pragma omp parallel for collapse(2)
   1432e:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   14332:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   14336:	8b 45 c4             	mov    -0x3c(%rbp),%eax
   14339:	89 45 f0             	mov    %eax,-0x10(%rbp)
   1433c:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14340:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   14344:	48 8d 45 e0          	lea    -0x20(%rbp),%rax
   14348:	b9 00 00 00 00       	mov    $0x0,%ecx
   1434d:	ba 00 00 00 00       	mov    $0x0,%edx
   14352:	48 89 c6             	mov    %rax,%rsi
   14355:	48 8d 05 a9 14 00 00 	lea    0x14a9(%rip),%rax        # 15805 <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.0>
   1435c:	48 89 c7             	mov    %rax,%rdi
   1435f:	e8 4c 64 ff ff       	call   a7b0 <GOMP_parallel@plt>
        for(int nn=0; nn<localParticlesX->ndim; nn++)
            *((localParticlesX->acceleration) + n*localParticlesX->ndim + nn) = 
                    *(recvBuf1 + n*localParticlesX->ndim + nn);
    
    double* recvBuf2;
    recvBuf2 = static_cast<double*>(malloc(localNParticlesY*(localParticlesY->ndim)*sizeof(double)));
   14364:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14368:	8b 40 30             	mov    0x30(%rax),%eax
   1436b:	0f af 45 c0          	imul   -0x40(%rbp),%eax
   1436f:	48 98                	cltq   
   14371:	48 c1 e0 03          	shl    $0x3,%rax
   14375:	48 89 c7             	mov    %rax,%rdi
   14378:	e8 b3 62 ff ff       	call   a630 <malloc@plt>
   1437d:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
    MPI_Allreduce(
   14381:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   14385:	48 8b 48 10          	mov    0x10(%rax),%rcx
        localParticlesY->acceleration,
        recvBuf2,
        localNParticlesY*localParticlesY->ndim,
   14389:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   1438d:	8b 40 30             	mov    0x30(%rax),%eax
    MPI_Allreduce(
   14390:	0f af 45 c0          	imul   -0x40(%rbp),%eax
   14394:	89 c2                	mov    %eax,%edx
        localParticlesY->acceleration,
   14396:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   1439a:	48 8b 40 18          	mov    0x18(%rax),%rax
    MPI_Allreduce(
   1439e:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
   143a2:	49 89 c9             	mov    %rcx,%r9
   143a5:	4c 8d 05 74 8c 00 00 	lea    0x8c74(%rip),%r8        # 1d020 <ompi_mpi_op_sum@@Base>
   143ac:	48 8d 0d 6d 98 00 00 	lea    0x986d(%rip),%rcx        # 1dc20 <ompi_mpi_double@@Base>
   143b3:	48 89 c7             	mov    %rax,%rdi
   143b6:	e8 a5 61 ff ff       	call   a560 <MPI_Allreduce@plt>
        MPI_DOUBLE,
        MPI_SUM,
        topology->reduceXComm
    );
    #pragma omp parallel for collapse(2)
   143bb:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   143bf:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   143c3:	8b 45 c0             	mov    -0x40(%rbp),%eax
   143c6:	89 45 f0             	mov    %eax,-0x10(%rbp)
   143c9:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   143cd:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   143d1:	48 8d 45 e0          	lea    -0x20(%rbp),%rax
   143d5:	b9 00 00 00 00       	mov    $0x0,%ecx
   143da:	ba 00 00 00 00       	mov    $0x0,%edx
   143df:	48 89 c6             	mov    %rax,%rsi
   143e2:	48 8d 05 61 15 00 00 	lea    0x1561(%rip),%rax        # 1594a <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.1>
   143e9:	48 89 c7             	mov    %rax,%rdi
   143ec:	e8 bf 63 ff ff       	call   a7b0 <GOMP_parallel@plt>
    for(int n=0; n<localNParticlesY; n++)
        for(int nn=0; nn<localParticlesY->ndim; nn++)
            *((localParticlesY->acceleration)+n*localParticlesY->ndim+nn) = 
                    *(recvBuf2+n*localParticlesY->ndim+nn);
    free(recvBuf1);
   143f1:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   143f5:	48 89 c7             	mov    %rax,%rdi
   143f8:	e8 03 61 ff ff       	call   a500 <free@plt>
    free(recvBuf2);
   143fd:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   14401:	48 89 c7             	mov    %rax,%rdi
   14404:	e8 f7 60 ff ff       	call   a500 <free@plt>
}
   14409:	90                   	nop
   1440a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1440e:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   14415:	00 00 
   14417:	74 05                	je     1441e <_Z6reduceP10particle_tiS0_iP10topology_t+0x174>
   14419:	e8 c2 60 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   1441e:	c9                   	leave  
   1441f:	c3                   	ret    

0000000000014420 <_Z6reduceP17chunk_particles_tS0_P10topology_t>:

void reduce(chunk_particles_t* localChunkParticlesX, chunk_particles_t* localChunkParticlesY, topology_t* topology) {
   14420:	f3 0f 1e fa          	endbr64 
   14424:	55                   	push   %rbp
   14425:	48 89 e5             	mov    %rsp,%rbp
   14428:	48 83 ec 20          	sub    $0x20,%rsp
   1442c:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   14430:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   14434:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    reduce(localChunkParticlesX->particles, localChunkParticlesX->nParticle, 
   14438:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   1443c:	8b 48 08             	mov    0x8(%rax),%ecx
   1443f:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   14443:	48 8b 10             	mov    (%rax),%rdx
   14446:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1444a:	8b 70 08             	mov    0x8(%rax),%esi
   1444d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14451:	48 8b 00             	mov    (%rax),%rax
   14454:	48 8b 7d e8          	mov    -0x18(%rbp),%rdi
   14458:	49 89 f8             	mov    %rdi,%r8
   1445b:	48 89 c7             	mov    %rax,%rdi
   1445e:	e8 47 fe ff ff       	call   142aa <_Z6reduceP10particle_tiS0_iP10topology_t>
            localChunkParticlesY->particles, localChunkParticlesY->nParticle, topology);
}
   14463:	90                   	nop
   14464:	c9                   	leave  
   14465:	c3                   	ret    

0000000000014466 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t>:


void gather(particle_t* localParticlesX, int localNParticlesX, 
                particle_t* localParticlesY, int localNParticlesY,
                particle_t** globalParticles, int* globalNParticles, topology_t* topology) {
   14466:	f3 0f 1e fa          	endbr64 
   1446a:	55                   	push   %rbp
   1446b:	48 89 e5             	mov    %rsp,%rbp
   1446e:	48 81 ec b0 00 00 00 	sub    $0xb0,%rsp
   14475:	48 89 bd 78 ff ff ff 	mov    %rdi,-0x88(%rbp)
   1447c:	89 b5 74 ff ff ff    	mov    %esi,-0x8c(%rbp)
   14482:	48 89 95 68 ff ff ff 	mov    %rdx,-0x98(%rbp)
   14489:	89 8d 70 ff ff ff    	mov    %ecx,-0x90(%rbp)
   1448f:	4c 89 85 60 ff ff ff 	mov    %r8,-0xa0(%rbp)
   14496:	4c 89 8d 58 ff ff ff 	mov    %r9,-0xa8(%rbp)
   1449d:	48 8b 45 10          	mov    0x10(%rbp),%rax
   144a1:	48 89 85 50 ff ff ff 	mov    %rax,-0xb0(%rbp)
   144a8:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   144af:	00 00 
   144b1:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   144b5:	31 c0                	xor    %eax,%eax
    // bool wholeParticlesInitNull = false;
    // if(*globalParticles==nullptr) {
        // wholeParticlesInitNull = true;
        // alloc_particles(globalParticles, 1, 1, 1);
    // }
    const int& rankX = topology->rankX;
   144b7:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   144be:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
    const int& rankY = topology->rankY;
   144c2:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   144c9:	48 83 c0 04          	add    $0x4,%rax
   144cd:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
    const int& gridX = topology->gridX;
   144d1:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   144d8:	48 83 c0 08          	add    $0x8,%rax
   144dc:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
    const int& gridY = topology->gridY;
   144e0:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   144e7:	48 83 c0 0c          	add    $0xc,%rax
   144eb:	48 89 45 b8          	mov    %rax,-0x48(%rbp)

    if(rankX==0 && rankY==0) {
   144ef:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   144f3:	8b 00                	mov    (%rax),%eax
   144f5:	85 c0                	test   %eax,%eax
   144f7:	0f 85 92 00 00 00    	jne    1458f <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x129>
   144fd:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   14501:	8b 00                	mov    (%rax),%eax
   14503:	85 c0                	test   %eax,%eax
   14505:	0f 85 84 00 00 00    	jne    1458f <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x129>
        assert(*globalNParticles == localNParticlesX * gridX);
   1450b:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   14512:	8b 10                	mov    (%rax),%edx
   14514:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   14518:	8b 00                	mov    (%rax),%eax
   1451a:	0f af 85 74 ff ff ff 	imul   -0x8c(%rbp),%eax
   14521:	39 c2                	cmp    %eax,%edx
   14523:	74 28                	je     1454d <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0xe7>
   14525:	48 8d 05 ac 29 00 00 	lea    0x29ac(%rip),%rax        # 16ed8 <_ZN3MPIL11LOCK_SHAREDE+0x238>
   1452c:	48 89 c1             	mov    %rax,%rcx
   1452f:	ba 3c 02 00 00       	mov    $0x23c,%edx
   14534:	48 8d 05 9d 27 00 00 	lea    0x279d(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   1453b:	48 89 c6             	mov    %rax,%rsi
   1453e:	48 8d 05 eb 29 00 00 	lea    0x29eb(%rip),%rax        # 16f30 <_ZN3MPIL11LOCK_SHAREDE+0x290>
   14545:	48 89 c7             	mov    %rax,%rdi
   14548:	e8 03 5c ff ff       	call   a150 <__assert_fail@plt>
        assert(*globalNParticles == localNParticlesY * gridY);
   1454d:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   14554:	8b 10                	mov    (%rax),%edx
   14556:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   1455a:	8b 00                	mov    (%rax),%eax
   1455c:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   14563:	39 c2                	cmp    %eax,%edx
   14565:	74 28                	je     1458f <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x129>
   14567:	48 8d 05 6a 29 00 00 	lea    0x296a(%rip),%rax        # 16ed8 <_ZN3MPIL11LOCK_SHAREDE+0x238>
   1456e:	48 89 c1             	mov    %rax,%rcx
   14571:	ba 3d 02 00 00       	mov    $0x23d,%edx
   14576:	48 8d 05 5b 27 00 00 	lea    0x275b(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   1457d:	48 89 c6             	mov    %rax,%rsi
   14580:	48 8d 05 d9 29 00 00 	lea    0x29d9(%rip),%rax        # 16f60 <_ZN3MPIL11LOCK_SHAREDE+0x2c0>
   14587:	48 89 c7             	mov    %rax,%rdi
   1458a:	e8 c1 5b ff ff       	call   a150 <__assert_fail@plt>
    }

    double *pos, *vel, *acc, *feat;
    if(rankX==0) {
   1458f:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   14593:	8b 00                	mov    (%rax),%eax
   14595:	85 c0                	test   %eax,%eax
   14597:	0f 85 2a 04 00 00    	jne    149c7 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x561>
        if(rankY==0) {
   1459d:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   145a1:	8b 00                	mov    (%rax),%eax
   145a3:	85 c0                	test   %eax,%eax
   145a5:	0f 85 91 02 00 00    	jne    1483c <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x3d6>
            pos = (*globalParticles)->position;
   145ab:	48 8b 85 60 ff ff ff 	mov    -0xa0(%rbp),%rax
   145b2:	48 8b 00             	mov    (%rax),%rax
   145b5:	48 8b 40 08          	mov    0x8(%rax),%rax
   145b9:	48 89 45 80          	mov    %rax,-0x80(%rbp)
            vel = (*globalParticles)->velocity;
   145bd:	48 8b 85 60 ff ff ff 	mov    -0xa0(%rbp),%rax
   145c4:	48 8b 00             	mov    (%rax),%rax
   145c7:	48 8b 40 10          	mov    0x10(%rax),%rax
   145cb:	48 89 45 88          	mov    %rax,-0x78(%rbp)
            acc = (*globalParticles)->acceleration;
   145cf:	48 8b 85 60 ff ff ff 	mov    -0xa0(%rbp),%rax
   145d6:	48 8b 00             	mov    (%rax),%rax
   145d9:	48 8b 40 18          	mov    0x18(%rax),%rax
   145dd:	48 89 45 90          	mov    %rax,-0x70(%rbp)
            feat = (*globalParticles)->features;
   145e1:	48 8b 85 60 ff ff ff 	mov    -0xa0(%rbp),%rax
   145e8:	48 8b 00             	mov    (%rax),%rax
   145eb:	48 8b 40 20          	mov    0x20(%rax),%rax
   145ef:	48 89 45 98          	mov    %rax,-0x68(%rbp)

            MPI_Gather(localParticlesY, 
   145f3:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   145fa:	48 8b 48 28          	mov    0x28(%rax),%rcx
                localNParticlesY*sizeof(particle_t),
                MPI_BYTE,
                *globalParticles,
                *globalNParticles*sizeof(particle_t)/gridY,
   145fe:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   14605:	8b 00                	mov    (%rax),%eax
   14607:	48 98                	cltq   
   14609:	48 c1 e0 06          	shl    $0x6,%rax
   1460d:	48 89 c6             	mov    %rax,%rsi
   14610:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14614:	8b 00                	mov    (%rax),%eax
   14616:	48 63 f8             	movslq %eax,%rdi
   14619:	48 89 f0             	mov    %rsi,%rax
   1461c:	ba 00 00 00 00       	mov    $0x0,%edx
   14621:	48 f7 f7             	div    %rdi
            MPI_Gather(localParticlesY, 
   14624:	89 c7                	mov    %eax,%edi
   14626:	48 8b 85 60 ff ff ff 	mov    -0xa0(%rbp),%rax
   1462d:	48 8b 10             	mov    (%rax),%rdx
                localNParticlesY*sizeof(particle_t),
   14630:	8b 85 70 ff ff ff    	mov    -0x90(%rbp),%eax
   14636:	48 98                	cltq   
            MPI_Gather(localParticlesY, 
   14638:	c1 e0 06             	shl    $0x6,%eax
   1463b:	89 c6                	mov    %eax,%esi
   1463d:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   14644:	51                   	push   %rcx
   14645:	6a 00                	push   $0x0
   14647:	4c 8d 0d d2 91 00 00 	lea    0x91d2(%rip),%r9        # 1d820 <ompi_mpi_byte@@Base>
   1464e:	41 89 f8             	mov    %edi,%r8d
   14651:	48 89 d1             	mov    %rdx,%rcx
   14654:	48 8d 15 c5 91 00 00 	lea    0x91c5(%rip),%rdx        # 1d820 <ompi_mpi_byte@@Base>
   1465b:	48 89 c7             	mov    %rax,%rdi
   1465e:	e8 2d 5c ff ff       	call   a290 <MPI_Gather@plt>
   14663:	48 83 c4 10          	add    $0x10,%rsp
                MPI_BYTE,
                0, 
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->position,
   14667:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   1466e:	48 8b 48 28          	mov    0x28(%rax),%rcx
                localNParticlesY*localParticlesY->ndim,
                MPI_DOUBLE,
                pos,
                (*globalNParticles)*(*globalParticles)->ndim/gridY,
   14672:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   14679:	8b 10                	mov    (%rax),%edx
   1467b:	48 8b 85 60 ff ff ff 	mov    -0xa0(%rbp),%rax
   14682:	48 8b 00             	mov    (%rax),%rax
   14685:	8b 40 30             	mov    0x30(%rax),%eax
   14688:	0f af c2             	imul   %edx,%eax
   1468b:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
   1468f:	8b 3a                	mov    (%rdx),%edi
            MPI_Gather(localParticlesY->position,
   14691:	99                   	cltd   
   14692:	f7 ff                	idiv   %edi
   14694:	89 c7                	mov    %eax,%edi
                localNParticlesY*localParticlesY->ndim,
   14696:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   1469d:	8b 40 30             	mov    0x30(%rax),%eax
            MPI_Gather(localParticlesY->position,
   146a0:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   146a7:	89 c6                	mov    %eax,%esi
   146a9:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   146b0:	48 8b 40 08          	mov    0x8(%rax),%rax
   146b4:	48 8b 55 80          	mov    -0x80(%rbp),%rdx
   146b8:	51                   	push   %rcx
   146b9:	6a 00                	push   $0x0
   146bb:	4c 8d 0d 5e 95 00 00 	lea    0x955e(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   146c2:	41 89 f8             	mov    %edi,%r8d
   146c5:	48 89 d1             	mov    %rdx,%rcx
   146c8:	48 8d 15 51 95 00 00 	lea    0x9551(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   146cf:	48 89 c7             	mov    %rax,%rdi
   146d2:	e8 b9 5b ff ff       	call   a290 <MPI_Gather@plt>
   146d7:	48 83 c4 10          	add    $0x10,%rsp
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->velocity,
   146db:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   146e2:	48 8b 48 28          	mov    0x28(%rax),%rcx
                localNParticlesY*localParticlesY->ndim,
                MPI_DOUBLE,
                vel,
                (*globalNParticles)*(*globalParticles)->ndim/gridY,
   146e6:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   146ed:	8b 10                	mov    (%rax),%edx
   146ef:	48 8b 85 60 ff ff ff 	mov    -0xa0(%rbp),%rax
   146f6:	48 8b 00             	mov    (%rax),%rax
   146f9:	8b 40 30             	mov    0x30(%rax),%eax
   146fc:	0f af c2             	imul   %edx,%eax
   146ff:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
   14703:	8b 3a                	mov    (%rdx),%edi
            MPI_Gather(localParticlesY->velocity,
   14705:	99                   	cltd   
   14706:	f7 ff                	idiv   %edi
   14708:	89 c7                	mov    %eax,%edi
                localNParticlesY*localParticlesY->ndim,
   1470a:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   14711:	8b 40 30             	mov    0x30(%rax),%eax
            MPI_Gather(localParticlesY->velocity,
   14714:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   1471b:	89 c6                	mov    %eax,%esi
   1471d:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   14724:	48 8b 40 10          	mov    0x10(%rax),%rax
   14728:	48 8b 55 88          	mov    -0x78(%rbp),%rdx
   1472c:	51                   	push   %rcx
   1472d:	6a 00                	push   $0x0
   1472f:	4c 8d 0d ea 94 00 00 	lea    0x94ea(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   14736:	41 89 f8             	mov    %edi,%r8d
   14739:	48 89 d1             	mov    %rdx,%rcx
   1473c:	48 8d 15 dd 94 00 00 	lea    0x94dd(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   14743:	48 89 c7             	mov    %rax,%rdi
   14746:	e8 45 5b ff ff       	call   a290 <MPI_Gather@plt>
   1474b:	48 83 c4 10          	add    $0x10,%rsp
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->acceleration,
   1474f:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   14756:	48 8b 48 28          	mov    0x28(%rax),%rcx
                localNParticlesY*localParticlesY->ndim,
                MPI_DOUBLE,
                acc,
                (*globalNParticles)*(*globalParticles)->ndim/gridY,
   1475a:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   14761:	8b 10                	mov    (%rax),%edx
   14763:	48 8b 85 60 ff ff ff 	mov    -0xa0(%rbp),%rax
   1476a:	48 8b 00             	mov    (%rax),%rax
   1476d:	8b 40 30             	mov    0x30(%rax),%eax
   14770:	0f af c2             	imul   %edx,%eax
   14773:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
   14777:	8b 3a                	mov    (%rdx),%edi
            MPI_Gather(localParticlesY->acceleration,
   14779:	99                   	cltd   
   1477a:	f7 ff                	idiv   %edi
   1477c:	89 c7                	mov    %eax,%edi
                localNParticlesY*localParticlesY->ndim,
   1477e:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   14785:	8b 40 30             	mov    0x30(%rax),%eax
            MPI_Gather(localParticlesY->acceleration,
   14788:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   1478f:	89 c6                	mov    %eax,%esi
   14791:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   14798:	48 8b 40 18          	mov    0x18(%rax),%rax
   1479c:	48 8b 55 90          	mov    -0x70(%rbp),%rdx
   147a0:	51                   	push   %rcx
   147a1:	6a 00                	push   $0x0
   147a3:	4c 8d 0d 76 94 00 00 	lea    0x9476(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   147aa:	41 89 f8             	mov    %edi,%r8d
   147ad:	48 89 d1             	mov    %rdx,%rcx
   147b0:	48 8d 15 69 94 00 00 	lea    0x9469(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   147b7:	48 89 c7             	mov    %rax,%rdi
   147ba:	e8 d1 5a ff ff       	call   a290 <MPI_Gather@plt>
   147bf:	48 83 c4 10          	add    $0x10,%rsp
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->features,
   147c3:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   147ca:	48 8b 48 28          	mov    0x28(%rax),%rcx
                localNParticlesY*localParticlesY->nfeat,
                MPI_DOUBLE,
                feat,
                (*globalNParticles)*(*globalParticles)->nfeat/gridY,
   147ce:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   147d5:	8b 10                	mov    (%rax),%edx
   147d7:	48 8b 85 60 ff ff ff 	mov    -0xa0(%rbp),%rax
   147de:	48 8b 00             	mov    (%rax),%rax
   147e1:	8b 40 34             	mov    0x34(%rax),%eax
   147e4:	0f af c2             	imul   %edx,%eax
   147e7:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
   147eb:	8b 3a                	mov    (%rdx),%edi
            MPI_Gather(localParticlesY->features,
   147ed:	99                   	cltd   
   147ee:	f7 ff                	idiv   %edi
   147f0:	89 c7                	mov    %eax,%edi
                localNParticlesY*localParticlesY->nfeat,
   147f2:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   147f9:	8b 40 34             	mov    0x34(%rax),%eax
            MPI_Gather(localParticlesY->features,
   147fc:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   14803:	89 c6                	mov    %eax,%esi
   14805:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   1480c:	48 8b 40 20          	mov    0x20(%rax),%rax
   14810:	48 8b 55 98          	mov    -0x68(%rbp),%rdx
   14814:	51                   	push   %rcx
   14815:	6a 00                	push   $0x0
   14817:	4c 8d 0d 02 94 00 00 	lea    0x9402(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   1481e:	41 89 f8             	mov    %edi,%r8d
   14821:	48 89 d1             	mov    %rdx,%rcx
   14824:	48 8d 15 f5 93 00 00 	lea    0x93f5(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   1482b:	48 89 c7             	mov    %rax,%rdi
   1482e:	e8 5d 5a ff ff       	call   a290 <MPI_Gather@plt>
   14833:	48 83 c4 10          	add    $0x10,%rsp
   14837:	e9 8b 01 00 00       	jmp    149c7 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x561>
                0,
                topology->scatterYComm
            );
        }
        else {
            MPI_Gather(localParticlesY, 
   1483c:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   14843:	48 8b 50 28          	mov    0x28(%rax),%rdx
                localNParticlesY*sizeof(particle_t),
   14847:	8b 85 70 ff ff ff    	mov    -0x90(%rbp),%eax
   1484d:	48 98                	cltq   
            MPI_Gather(localParticlesY, 
   1484f:	c1 e0 06             	shl    $0x6,%eax
   14852:	89 c6                	mov    %eax,%esi
   14854:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   1485b:	52                   	push   %rdx
   1485c:	6a 00                	push   $0x0
   1485e:	4c 8d 0d bb 8f 00 00 	lea    0x8fbb(%rip),%r9        # 1d820 <ompi_mpi_byte@@Base>
   14865:	41 b8 ff ff ff ff    	mov    $0xffffffff,%r8d
   1486b:	b9 00 00 00 00       	mov    $0x0,%ecx
   14870:	48 8d 15 a9 8f 00 00 	lea    0x8fa9(%rip),%rdx        # 1d820 <ompi_mpi_byte@@Base>
   14877:	48 89 c7             	mov    %rax,%rdi
   1487a:	e8 11 5a ff ff       	call   a290 <MPI_Gather@plt>
   1487f:	48 83 c4 10          	add    $0x10,%rsp
                MPI_BYTE,
                0, 
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->position,
   14883:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   1488a:	48 8b 50 28          	mov    0x28(%rax),%rdx
                localNParticlesY*localParticlesY->ndim,
   1488e:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   14895:	8b 40 30             	mov    0x30(%rax),%eax
            MPI_Gather(localParticlesY->position,
   14898:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   1489f:	89 c6                	mov    %eax,%esi
   148a1:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   148a8:	48 8b 40 08          	mov    0x8(%rax),%rax
   148ac:	52                   	push   %rdx
   148ad:	6a 00                	push   $0x0
   148af:	4c 8d 0d 6a 93 00 00 	lea    0x936a(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   148b6:	41 b8 ff ff ff ff    	mov    $0xffffffff,%r8d
   148bc:	b9 00 00 00 00       	mov    $0x0,%ecx
   148c1:	48 8d 15 58 93 00 00 	lea    0x9358(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   148c8:	48 89 c7             	mov    %rax,%rdi
   148cb:	e8 c0 59 ff ff       	call   a290 <MPI_Gather@plt>
   148d0:	48 83 c4 10          	add    $0x10,%rsp
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->velocity,
   148d4:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   148db:	48 8b 50 28          	mov    0x28(%rax),%rdx
                localNParticlesY*localParticlesY->ndim,
   148df:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   148e6:	8b 40 30             	mov    0x30(%rax),%eax
            MPI_Gather(localParticlesY->velocity,
   148e9:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   148f0:	89 c6                	mov    %eax,%esi
   148f2:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   148f9:	48 8b 40 10          	mov    0x10(%rax),%rax
   148fd:	52                   	push   %rdx
   148fe:	6a 00                	push   $0x0
   14900:	4c 8d 0d 19 93 00 00 	lea    0x9319(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   14907:	41 b8 ff ff ff ff    	mov    $0xffffffff,%r8d
   1490d:	b9 00 00 00 00       	mov    $0x0,%ecx
   14912:	48 8d 15 07 93 00 00 	lea    0x9307(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   14919:	48 89 c7             	mov    %rax,%rdi
   1491c:	e8 6f 59 ff ff       	call   a290 <MPI_Gather@plt>
   14921:	48 83 c4 10          	add    $0x10,%rsp
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->acceleration,
   14925:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   1492c:	48 8b 50 28          	mov    0x28(%rax),%rdx
                localNParticlesY*localParticlesY->ndim,
   14930:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   14937:	8b 40 30             	mov    0x30(%rax),%eax
            MPI_Gather(localParticlesY->acceleration,
   1493a:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   14941:	89 c6                	mov    %eax,%esi
   14943:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   1494a:	48 8b 40 18          	mov    0x18(%rax),%rax
   1494e:	52                   	push   %rdx
   1494f:	6a 00                	push   $0x0
   14951:	4c 8d 0d c8 92 00 00 	lea    0x92c8(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   14958:	41 b8 ff ff ff ff    	mov    $0xffffffff,%r8d
   1495e:	b9 00 00 00 00       	mov    $0x0,%ecx
   14963:	48 8d 15 b6 92 00 00 	lea    0x92b6(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   1496a:	48 89 c7             	mov    %rax,%rdi
   1496d:	e8 1e 59 ff ff       	call   a290 <MPI_Gather@plt>
   14972:	48 83 c4 10          	add    $0x10,%rsp
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->features,
   14976:	48 8b 85 50 ff ff ff 	mov    -0xb0(%rbp),%rax
   1497d:	48 8b 50 28          	mov    0x28(%rax),%rdx
                localNParticlesY*localParticlesY->nfeat,
   14981:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   14988:	8b 40 34             	mov    0x34(%rax),%eax
            MPI_Gather(localParticlesY->features,
   1498b:	0f af 85 70 ff ff ff 	imul   -0x90(%rbp),%eax
   14992:	89 c6                	mov    %eax,%esi
   14994:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   1499b:	48 8b 40 20          	mov    0x20(%rax),%rax
   1499f:	52                   	push   %rdx
   149a0:	6a 00                	push   $0x0
   149a2:	4c 8d 0d 77 92 00 00 	lea    0x9277(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   149a9:	41 b8 ff ff ff ff    	mov    $0xffffffff,%r8d
   149af:	b9 00 00 00 00       	mov    $0x0,%ecx
   149b4:	48 8d 15 65 92 00 00 	lea    0x9265(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   149bb:	48 89 c7             	mov    %rax,%rdi
   149be:	e8 cd 58 ff ff       	call   a290 <MPI_Gather@plt>
   149c3:	48 83 c4 10          	add    $0x10,%rsp
                0,
                topology->scatterYComm
            );
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
   149c7:	48 8d 05 52 90 00 00 	lea    0x9052(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   149ce:	48 89 c7             	mov    %rax,%rdi
   149d1:	e8 0a 58 ff ff       	call   a1e0 <MPI_Barrier@plt>

    if(rankX==0 && rankY==0) {
   149d6:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   149da:	8b 00                	mov    (%rax),%eax
   149dc:	85 c0                	test   %eax,%eax
   149de:	75 61                	jne    14a41 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x5db>
   149e0:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   149e4:	8b 00                	mov    (%rax),%eax
   149e6:	85 c0                	test   %eax,%eax
   149e8:	75 57                	jne    14a41 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x5db>
        #pragma omp parallel for
   149ea:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   149ee:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   149f2:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   149f6:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   149fa:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   149fe:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   14a02:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   14a06:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   14a0a:	48 8b 85 58 ff ff ff 	mov    -0xa8(%rbp),%rax
   14a11:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   14a15:	48 8b 85 60 ff ff ff 	mov    -0xa0(%rbp),%rax
   14a1c:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   14a20:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   14a24:	b9 00 00 00 00       	mov    $0x0,%ecx
   14a29:	ba 00 00 00 00       	mov    $0x0,%edx
   14a2e:	48 89 c6             	mov    %rax,%rsi
   14a31:	48 8d 05 57 10 00 00 	lea    0x1057(%rip),%rax        # 15a8f <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t._omp_fn.0>
   14a38:	48 89 c7             	mov    %rax,%rdi
   14a3b:	e8 70 5d ff ff       	call   a7b0 <GOMP_parallel@plt>
            (*globalParticles)->velocity = vel+n*((*globalParticles)->ndim);
            (*globalParticles)->acceleration = acc+n*((*globalParticles)->ndim);
            (*globalParticles)->features = feat+n*((*globalParticles)->nfeat);
        }
    }
    return;
   14a40:	90                   	nop
   14a41:	90                   	nop
}
   14a42:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14a46:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   14a4d:	00 00 
   14a4f:	74 05                	je     14a56 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x5f0>
   14a51:	e8 8a 5a ff ff       	call   a4e0 <__stack_chk_fail@plt>
   14a56:	c9                   	leave  
   14a57:	c3                   	ret    

0000000000014a58 <_Z6gatherP17chunk_particles_tS0_S0_P10topology_t>:

void gather(chunk_particles_t* localChunkParticlesX, 
            chunk_particles_t* localChunkParticlesY, 
            chunk_particles_t* globalChunkParticles, 
            topology_t* topology) {
   14a58:	f3 0f 1e fa          	endbr64 
   14a5c:	55                   	push   %rbp
   14a5d:	48 89 e5             	mov    %rsp,%rbp
   14a60:	48 83 ec 50          	sub    $0x50,%rsp
   14a64:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
   14a68:	48 89 75 c0          	mov    %rsi,-0x40(%rbp)
   14a6c:	48 89 55 b8          	mov    %rdx,-0x48(%rbp)
   14a70:	48 89 4d b0          	mov    %rcx,-0x50(%rbp)
   14a74:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   14a7b:	00 00 
   14a7d:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   14a81:	31 c0                	xor    %eax,%eax
    particle_t* localParticlesX = localChunkParticlesX->particles;
   14a83:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14a87:	48 8b 00             	mov    (%rax),%rax
   14a8a:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    int localNParticlesX = localChunkParticlesX->nParticle;
   14a8e:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14a92:	8b 40 08             	mov    0x8(%rax),%eax
   14a95:	89 45 d8             	mov    %eax,-0x28(%rbp)
    particle_t* localParticlesY = localChunkParticlesY->particles;
   14a98:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   14a9c:	48 8b 00             	mov    (%rax),%rax
   14a9f:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
    int localNParticlesY = localChunkParticlesY->nParticle;
   14aa3:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   14aa7:	8b 40 08             	mov    0x8(%rax),%eax
   14aaa:	89 45 dc             	mov    %eax,-0x24(%rbp)
    particle_t* globalParticles;
    int globalNParticles;
    if(globalChunkParticles==nullptr) {
   14aad:	48 83 7d b8 00       	cmpq   $0x0,-0x48(%rbp)
   14ab2:	75 11                	jne    14ac5 <_Z6gatherP17chunk_particles_tS0_S0_P10topology_t+0x6d>
        globalParticles = nullptr;
   14ab4:	48 c7 45 e0 00 00 00 	movq   $0x0,-0x20(%rbp)
   14abb:	00 
        globalNParticles = -1;
   14abc:	c7 45 d4 ff ff ff ff 	movl   $0xffffffff,-0x2c(%rbp)
   14ac3:	eb 45                	jmp    14b0a <_Z6gatherP17chunk_particles_tS0_S0_P10topology_t+0xb2>
    }
    else {
        globalParticles = globalChunkParticles->particles;
   14ac5:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14ac9:	48 8b 00             	mov    (%rax),%rax
   14acc:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
        globalNParticles = globalChunkParticles->nParticle;
   14ad0:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14ad4:	8b 40 08             	mov    0x8(%rax),%eax
   14ad7:	89 45 d4             	mov    %eax,-0x2c(%rbp)
        printf("%s:%d %lx %d\n", __FILE__, __LINE__, globalParticles, globalNParticles);
   14ada:	8b 55 d4             	mov    -0x2c(%rbp),%edx
   14add:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   14ae1:	41 89 d0             	mov    %edx,%r8d
   14ae4:	48 89 c1             	mov    %rax,%rcx
   14ae7:	ba cd 02 00 00       	mov    $0x2cd,%edx
   14aec:	48 8d 05 e5 21 00 00 	lea    0x21e5(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   14af3:	48 89 c6             	mov    %rax,%rsi
   14af6:	48 8d 05 91 24 00 00 	lea    0x2491(%rip),%rax        # 16f8e <_ZN3MPIL11LOCK_SHAREDE+0x2ee>
   14afd:	48 89 c7             	mov    %rax,%rdi
   14b00:	b8 00 00 00 00       	mov    $0x0,%eax
   14b05:	e8 d6 52 ff ff       	call   9de0 <printf@plt>
    }

    gather(localParticlesX, localNParticlesX, localParticlesY, 
   14b0a:	4c 8d 45 d4          	lea    -0x2c(%rbp),%r8
   14b0e:	48 8d 7d e0          	lea    -0x20(%rbp),%rdi
   14b12:	8b 4d dc             	mov    -0x24(%rbp),%ecx
   14b15:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   14b19:	8b 75 d8             	mov    -0x28(%rbp),%esi
   14b1c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   14b20:	48 83 ec 08          	sub    $0x8,%rsp
   14b24:	ff 75 b0             	push   -0x50(%rbp)
   14b27:	4d 89 c1             	mov    %r8,%r9
   14b2a:	49 89 f8             	mov    %rdi,%r8
   14b2d:	48 89 c7             	mov    %rax,%rdi
   14b30:	e8 31 f9 ff ff       	call   14466 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t>
   14b35:	48 83 c4 10          	add    $0x10,%rsp
            localNParticlesY, &globalParticles, &globalNParticles, topology);
    if(globalChunkParticles!=nullptr) {
   14b39:	48 83 7d b8 00       	cmpq   $0x0,-0x48(%rbp)
   14b3e:	74 45                	je     14b85 <_Z6gatherP17chunk_particles_tS0_S0_P10topology_t+0x12d>
        printf("%s:%d %lx %d\n", __FILE__, __LINE__, globalParticles, globalNParticles);
   14b40:	8b 55 d4             	mov    -0x2c(%rbp),%edx
   14b43:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   14b47:	41 89 d0             	mov    %edx,%r8d
   14b4a:	48 89 c1             	mov    %rax,%rcx
   14b4d:	ba d3 02 00 00       	mov    $0x2d3,%edx
   14b52:	48 8d 05 7f 21 00 00 	lea    0x217f(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   14b59:	48 89 c6             	mov    %rax,%rsi
   14b5c:	48 8d 05 2b 24 00 00 	lea    0x242b(%rip),%rax        # 16f8e <_ZN3MPIL11LOCK_SHAREDE+0x2ee>
   14b63:	48 89 c7             	mov    %rax,%rdi
   14b66:	b8 00 00 00 00       	mov    $0x0,%eax
   14b6b:	e8 70 52 ff ff       	call   9de0 <printf@plt>
        globalChunkParticles->particles = globalParticles;
   14b70:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
   14b74:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14b78:	48 89 10             	mov    %rdx,(%rax)
        globalChunkParticles->nParticle = globalNParticles;
   14b7b:	8b 55 d4             	mov    -0x2c(%rbp),%edx
   14b7e:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14b82:	89 50 08             	mov    %edx,0x8(%rax)
    }
   14b85:	90                   	nop
   14b86:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14b8a:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   14b91:	00 00 
   14b93:	74 05                	je     14b9a <_Z6gatherP17chunk_particles_tS0_S0_P10topology_t+0x142>
   14b95:	e8 46 59 ff ff       	call   a4e0 <__stack_chk_fail@plt>
   14b9a:	c9                   	leave  
   14b9b:	c3                   	ret    

0000000000014b9c <_Z41__static_initialization_and_destruction_0ii>:
   14b9c:	f3 0f 1e fa          	endbr64 
   14ba0:	55                   	push   %rbp
   14ba1:	48 89 e5             	mov    %rsp,%rbp
   14ba4:	48 83 ec 10          	sub    $0x10,%rsp
   14ba8:	89 7d fc             	mov    %edi,-0x4(%rbp)
   14bab:	89 75 f8             	mov    %esi,-0x8(%rbp)
   14bae:	83 7d fc 01          	cmpl   $0x1,-0x4(%rbp)
   14bb2:	75 3b                	jne    14bef <_Z41__static_initialization_and_destruction_0ii+0x53>
   14bb4:	81 7d f8 ff ff 00 00 	cmpl   $0xffff,-0x8(%rbp)
   14bbb:	75 32                	jne    14bef <_Z41__static_initialization_and_destruction_0ii+0x53>
   14bbd:	48 8d 05 8f 97 00 00 	lea    0x978f(%rip),%rax        # 1e353 <_ZStL8__ioinit>
   14bc4:	48 89 c7             	mov    %rax,%rdi
   14bc7:	e8 04 5c ff ff       	call   a7d0 <_ZNSt8ios_base4InitC1Ev@plt>
   14bcc:	48 8d 05 35 84 00 00 	lea    0x8435(%rip),%rax        # 1d008 <__dso_handle>
   14bd3:	48 89 c2             	mov    %rax,%rdx
   14bd6:	48 8d 05 76 97 00 00 	lea    0x9776(%rip),%rax        # 1e353 <_ZStL8__ioinit>
   14bdd:	48 89 c6             	mov    %rax,%rsi
   14be0:	48 8b 05 11 84 00 00 	mov    0x8411(%rip),%rax        # 1cff8 <_ZNSt8ios_base4InitD1Ev@GLIBCXX_3.4>
   14be7:	48 89 c7             	mov    %rax,%rdi
   14bea:	e8 41 57 ff ff       	call   a330 <__cxa_atexit@plt>
   14bef:	90                   	nop
   14bf0:	c9                   	leave  
   14bf1:	c3                   	ret    

0000000000014bf2 <_GLOBAL__sub_I__Z18setupCommunicatorsP10topology_tii>:
   14bf2:	f3 0f 1e fa          	endbr64 
   14bf6:	55                   	push   %rbp
   14bf7:	48 89 e5             	mov    %rsp,%rbp
   14bfa:	be ff ff 00 00       	mov    $0xffff,%esi
   14bff:	bf 01 00 00 00       	mov    $0x1,%edi
   14c04:	e8 93 ff ff ff       	call   14b9c <_Z41__static_initialization_and_destruction_0ii>
   14c09:	5d                   	pop    %rbp
   14c0a:	c3                   	ret    

0000000000014c0b <_Z7paddingPP10particle_tPii._omp_fn.0>:
    #pragma omp parallel for collapse(2)
   14c0b:	f3 0f 1e fa          	endbr64 
   14c0f:	55                   	push   %rbp
   14c10:	48 89 e5             	mov    %rsp,%rbp
   14c13:	41 57                	push   %r15
   14c15:	41 56                	push   %r14
   14c17:	41 55                	push   %r13
   14c19:	41 54                	push   %r12
   14c1b:	53                   	push   %rbx
   14c1c:	48 83 ec 48          	sub    $0x48,%rsp
   14c20:	48 89 7d 98          	mov    %rdi,-0x68(%rbp)
   14c24:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   14c28:	8b 40 20             	mov    0x20(%rax),%eax
   14c2b:	89 45 ac             	mov    %eax,-0x54(%rbp)
   14c2e:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   14c32:	48 8b 40 10          	mov    0x10(%rax),%rax
   14c36:	48 89 45 b8          	mov    %rax,-0x48(%rbp)
   14c3a:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   14c3e:	48 8b 40 08          	mov    0x8(%rax),%rax
   14c42:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   14c46:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   14c4a:	48 8b 00             	mov    (%rax),%rax
   14c4d:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   14c51:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   14c55:	8b 40 1c             	mov    0x1c(%rax),%eax
   14c58:	89 45 b0             	mov    %eax,-0x50(%rbp)
   14c5b:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   14c5f:	8b 40 18             	mov    0x18(%rax),%eax
   14c62:	89 45 b4             	mov    %eax,-0x4c(%rbp)
    for(int n=oriSize; n<newSize; n++) 
   14c65:	44 8b 75 b4          	mov    -0x4c(%rbp),%r14d
   14c69:	8b 45 b0             	mov    -0x50(%rbp),%eax
   14c6c:	8b 5d ac             	mov    -0x54(%rbp),%ebx
   14c6f:	41 39 c6             	cmp    %eax,%r14d
   14c72:	0f 8d 23 01 00 00    	jge    14d9b <_Z7paddingPP10particle_tPii._omp_fn.0+0x190>
   14c78:	44 29 f0             	sub    %r14d,%eax
   14c7b:	85 db                	test   %ebx,%ebx
   14c7d:	0f 8e 18 01 00 00    	jle    14d9b <_Z7paddingPP10particle_tPii._omp_fn.0+0x190>
   14c83:	41 89 dd             	mov    %ebx,%r13d
   14c86:	41 0f af c5          	imul   %r13d,%eax
   14c8a:	41 89 c4             	mov    %eax,%r12d
   14c8d:	e8 ae 5d ff ff       	call   aa40 <omp_get_num_threads@plt>
   14c92:	41 89 c7             	mov    %eax,%r15d
   14c95:	e8 e6 5b ff ff       	call   a880 <omp_get_thread_num@plt>
   14c9a:	89 c1                	mov    %eax,%ecx
   14c9c:	44 89 e0             	mov    %r12d,%eax
   14c9f:	ba 00 00 00 00       	mov    $0x0,%edx
   14ca4:	41 f7 f7             	div    %r15d
   14ca7:	89 c6                	mov    %eax,%esi
   14ca9:	44 89 e0             	mov    %r12d,%eax
   14cac:	ba 00 00 00 00       	mov    $0x0,%edx
   14cb1:	41 f7 f7             	div    %r15d
   14cb4:	89 d0                	mov    %edx,%eax
   14cb6:	39 c1                	cmp    %eax,%ecx
   14cb8:	0f 82 d0 00 00 00    	jb     14d8e <_Z7paddingPP10particle_tPii._omp_fn.0+0x183>
   14cbe:	0f af ce             	imul   %esi,%ecx
   14cc1:	89 ca                	mov    %ecx,%edx
   14cc3:	8d 0c 02             	lea    (%rdx,%rax,1),%ecx
   14cc6:	01 ce                	add    %ecx,%esi
   14cc8:	39 f1                	cmp    %esi,%ecx
   14cca:	0f 83 cb 00 00 00    	jae    14d9b <_Z7paddingPP10particle_tPii._omp_fn.0+0x190>
   14cd0:	89 c8                	mov    %ecx,%eax
   14cd2:	ba 00 00 00 00       	mov    $0x0,%edx
   14cd7:	41 f7 f5             	div    %r13d
   14cda:	89 d0                	mov    %edx,%eax
   14cdc:	89 45 a4             	mov    %eax,-0x5c(%rbp)
   14cdf:	89 c8                	mov    %ecx,%eax
   14ce1:	ba 00 00 00 00       	mov    $0x0,%edx
   14ce6:	41 f7 f5             	div    %r13d
   14ce9:	44 01 f0             	add    %r14d,%eax
   14cec:	89 45 a8             	mov    %eax,-0x58(%rbp)
            *(newPos+n*ndim+nn) = 0;
   14cef:	8b 45 a8             	mov    -0x58(%rbp),%eax
   14cf2:	0f af 45 ac          	imul   -0x54(%rbp),%eax
   14cf6:	48 63 d0             	movslq %eax,%rdx
   14cf9:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   14cfc:	48 98                	cltq   
   14cfe:	48 01 d0             	add    %rdx,%rax
   14d01:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   14d08:	00 
   14d09:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14d0d:	48 01 d0             	add    %rdx,%rax
   14d10:	66 0f ef c0          	pxor   %xmm0,%xmm0
   14d14:	f2 0f 11 00          	movsd  %xmm0,(%rax)
            *(newVel+n*ndim+nn) = 0;
   14d18:	8b 45 a8             	mov    -0x58(%rbp),%eax
   14d1b:	0f af 45 ac          	imul   -0x54(%rbp),%eax
   14d1f:	48 63 d0             	movslq %eax,%rdx
   14d22:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   14d25:	48 98                	cltq   
   14d27:	48 01 d0             	add    %rdx,%rax
   14d2a:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   14d31:	00 
   14d32:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   14d36:	48 01 d0             	add    %rdx,%rax
   14d39:	66 0f ef c0          	pxor   %xmm0,%xmm0
   14d3d:	f2 0f 11 00          	movsd  %xmm0,(%rax)
            *(newAcc+n*ndim+nn) = 0;
   14d41:	8b 45 a8             	mov    -0x58(%rbp),%eax
   14d44:	0f af 45 ac          	imul   -0x54(%rbp),%eax
   14d48:	48 63 d0             	movslq %eax,%rdx
   14d4b:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   14d4e:	48 98                	cltq   
   14d50:	48 01 d0             	add    %rdx,%rax
   14d53:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   14d5a:	00 
   14d5b:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14d5f:	48 01 d0             	add    %rdx,%rax
   14d62:	66 0f ef c0          	pxor   %xmm0,%xmm0
   14d66:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   14d6a:	83 c1 01             	add    $0x1,%ecx
   14d6d:	39 f1                	cmp    %esi,%ecx
   14d6f:	73 2a                	jae    14d9b <_Z7paddingPP10particle_tPii._omp_fn.0+0x190>
   14d71:	83 45 a4 01          	addl   $0x1,-0x5c(%rbp)
   14d75:	39 5d a4             	cmp    %ebx,-0x5c(%rbp)
   14d78:	0f 8c 71 ff ff ff    	jl     14cef <_Z7paddingPP10particle_tPii._omp_fn.0+0xe4>
   14d7e:	c7 45 a4 00 00 00 00 	movl   $0x0,-0x5c(%rbp)
   14d85:	83 45 a8 01          	addl   $0x1,-0x58(%rbp)
   14d89:	e9 61 ff ff ff       	jmp    14cef <_Z7paddingPP10particle_tPii._omp_fn.0+0xe4>
   14d8e:	b8 00 00 00 00       	mov    $0x0,%eax
   14d93:	83 c6 01             	add    $0x1,%esi
   14d96:	e9 23 ff ff ff       	jmp    14cbe <_Z7paddingPP10particle_tPii._omp_fn.0+0xb3>
    #pragma omp parallel for collapse(2)
   14d9b:	90                   	nop
   14d9c:	48 83 c4 48          	add    $0x48,%rsp
   14da0:	5b                   	pop    %rbx
   14da1:	41 5c                	pop    %r12
   14da3:	41 5d                	pop    %r13
   14da5:	41 5e                	pop    %r14
   14da7:	41 5f                	pop    %r15
   14da9:	5d                   	pop    %rbp
   14daa:	c3                   	ret    

0000000000014dab <_Z7paddingPP10particle_tPii._omp_fn.1>:
    #pragma omp parallel for collapse(2)
   14dab:	f3 0f 1e fa          	endbr64 
   14daf:	55                   	push   %rbp
   14db0:	48 89 e5             	mov    %rsp,%rbp
   14db3:	41 57                	push   %r15
   14db5:	41 56                	push   %r14
   14db7:	41 55                	push   %r13
   14db9:	41 54                	push   %r12
   14dbb:	53                   	push   %rbx
   14dbc:	48 83 ec 38          	sub    $0x38,%rsp
   14dc0:	48 89 7d a8          	mov    %rdi,-0x58(%rbp)
   14dc4:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   14dc8:	8b 40 10             	mov    0x10(%rax),%eax
   14dcb:	89 45 bc             	mov    %eax,-0x44(%rbp)
   14dce:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   14dd2:	48 8b 00             	mov    (%rax),%rax
   14dd5:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   14dd9:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   14ddd:	8b 40 0c             	mov    0xc(%rax),%eax
   14de0:	89 45 c0             	mov    %eax,-0x40(%rbp)
   14de3:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   14de7:	8b 40 08             	mov    0x8(%rax),%eax
   14dea:	89 45 c4             	mov    %eax,-0x3c(%rbp)
    for(int n=oriSize; n<newSize; n++)
   14ded:	44 8b 75 c4          	mov    -0x3c(%rbp),%r14d
   14df1:	8b 45 c0             	mov    -0x40(%rbp),%eax
   14df4:	8b 5d bc             	mov    -0x44(%rbp),%ebx
   14df7:	41 39 c6             	cmp    %eax,%r14d
   14dfa:	0f 8d bf 00 00 00    	jge    14ebf <_Z7paddingPP10particle_tPii._omp_fn.1+0x114>
   14e00:	44 29 f0             	sub    %r14d,%eax
   14e03:	85 db                	test   %ebx,%ebx
   14e05:	0f 8e b4 00 00 00    	jle    14ebf <_Z7paddingPP10particle_tPii._omp_fn.1+0x114>
   14e0b:	41 89 dd             	mov    %ebx,%r13d
   14e0e:	41 0f af c5          	imul   %r13d,%eax
   14e12:	41 89 c4             	mov    %eax,%r12d
   14e15:	e8 26 5c ff ff       	call   aa40 <omp_get_num_threads@plt>
   14e1a:	41 89 c7             	mov    %eax,%r15d
   14e1d:	e8 5e 5a ff ff       	call   a880 <omp_get_thread_num@plt>
   14e22:	89 c1                	mov    %eax,%ecx
   14e24:	44 89 e0             	mov    %r12d,%eax
   14e27:	ba 00 00 00 00       	mov    $0x0,%edx
   14e2c:	41 f7 f7             	div    %r15d
   14e2f:	89 c6                	mov    %eax,%esi
   14e31:	44 89 e0             	mov    %r12d,%eax
   14e34:	ba 00 00 00 00       	mov    $0x0,%edx
   14e39:	41 f7 f7             	div    %r15d
   14e3c:	89 d0                	mov    %edx,%eax
   14e3e:	39 c1                	cmp    %eax,%ecx
   14e40:	72 73                	jb     14eb5 <_Z7paddingPP10particle_tPii._omp_fn.1+0x10a>
   14e42:	0f af ce             	imul   %esi,%ecx
   14e45:	89 ca                	mov    %ecx,%edx
   14e47:	8d 0c 02             	lea    (%rdx,%rax,1),%ecx
   14e4a:	01 ce                	add    %ecx,%esi
   14e4c:	39 f1                	cmp    %esi,%ecx
   14e4e:	73 6f                	jae    14ebf <_Z7paddingPP10particle_tPii._omp_fn.1+0x114>
   14e50:	89 c8                	mov    %ecx,%eax
   14e52:	ba 00 00 00 00       	mov    $0x0,%edx
   14e57:	41 f7 f5             	div    %r13d
   14e5a:	89 d0                	mov    %edx,%eax
   14e5c:	89 45 b4             	mov    %eax,-0x4c(%rbp)
   14e5f:	89 c8                	mov    %ecx,%eax
   14e61:	ba 00 00 00 00       	mov    $0x0,%edx
   14e66:	41 f7 f5             	div    %r13d
   14e69:	44 01 f0             	add    %r14d,%eax
   14e6c:	89 45 b8             	mov    %eax,-0x48(%rbp)
            *(newFeat+n*nfeat+nn) = 0;
   14e6f:	8b 45 b8             	mov    -0x48(%rbp),%eax
   14e72:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   14e76:	48 63 d0             	movslq %eax,%rdx
   14e79:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   14e7c:	48 98                	cltq   
   14e7e:	48 01 d0             	add    %rdx,%rax
   14e81:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   14e88:	00 
   14e89:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14e8d:	48 01 d0             	add    %rdx,%rax
   14e90:	66 0f ef c0          	pxor   %xmm0,%xmm0
   14e94:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   14e98:	83 c1 01             	add    $0x1,%ecx
   14e9b:	39 f1                	cmp    %esi,%ecx
   14e9d:	73 20                	jae    14ebf <_Z7paddingPP10particle_tPii._omp_fn.1+0x114>
   14e9f:	83 45 b4 01          	addl   $0x1,-0x4c(%rbp)
   14ea3:	39 5d b4             	cmp    %ebx,-0x4c(%rbp)
   14ea6:	7c c7                	jl     14e6f <_Z7paddingPP10particle_tPii._omp_fn.1+0xc4>
   14ea8:	c7 45 b4 00 00 00 00 	movl   $0x0,-0x4c(%rbp)
   14eaf:	83 45 b8 01          	addl   $0x1,-0x48(%rbp)
   14eb3:	eb ba                	jmp    14e6f <_Z7paddingPP10particle_tPii._omp_fn.1+0xc4>
   14eb5:	b8 00 00 00 00       	mov    $0x0,%eax
   14eba:	83 c6 01             	add    $0x1,%esi
   14ebd:	eb 83                	jmp    14e42 <_Z7paddingPP10particle_tPii._omp_fn.1+0x97>
    #pragma omp parallel for collapse(2)
   14ebf:	90                   	nop
   14ec0:	48 83 c4 38          	add    $0x38,%rsp
   14ec4:	5b                   	pop    %rbx
   14ec5:	41 5c                	pop    %r12
   14ec7:	41 5d                	pop    %r13
   14ec9:	41 5e                	pop    %r14
   14ecb:	41 5f                	pop    %r15
   14ecd:	5d                   	pop    %rbp
   14ece:	c3                   	ret    

0000000000014ecf <_Z7paddingPP10particle_tPii._omp_fn.2>:
    #pragma omp parallel for
   14ecf:	f3 0f 1e fa          	endbr64 
   14ed3:	55                   	push   %rbp
   14ed4:	48 89 e5             	mov    %rsp,%rbp
   14ed7:	41 54                	push   %r12
   14ed9:	53                   	push   %rbx
   14eda:	48 83 ec 70          	sub    $0x70,%rsp
   14ede:	48 89 7d 88          	mov    %rdi,-0x78(%rbp)
   14ee2:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14ee6:	8b 40 54             	mov    0x54(%rax),%eax
   14ee9:	89 45 98             	mov    %eax,-0x68(%rbp)
   14eec:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14ef0:	8b 40 50             	mov    0x50(%rax),%eax
   14ef3:	89 45 9c             	mov    %eax,-0x64(%rbp)
   14ef6:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14efa:	48 8b 40 40          	mov    0x40(%rax),%rax
   14efe:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
   14f02:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14f06:	48 8b 40 38          	mov    0x38(%rax),%rax
   14f0a:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
   14f0e:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14f12:	48 8b 40 30          	mov    0x30(%rax),%rax
   14f16:	48 89 45 b8          	mov    %rax,-0x48(%rbp)
   14f1a:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14f1e:	48 8b 40 28          	mov    0x28(%rax),%rax
   14f22:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   14f26:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14f2a:	48 8b 40 20          	mov    0x20(%rax),%rax
   14f2e:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   14f32:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14f36:	48 8b 40 18          	mov    0x18(%rax),%rax
   14f3a:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   14f3e:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14f42:	48 8b 40 10          	mov    0x10(%rax),%rax
   14f46:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   14f4a:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14f4e:	48 8b 40 08          	mov    0x8(%rax),%rax
   14f52:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   14f56:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14f5a:	48 8b 00             	mov    (%rax),%rax
   14f5d:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   14f61:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14f65:	8b 40 4c             	mov    0x4c(%rax),%eax
   14f68:	89 45 a0             	mov    %eax,-0x60(%rbp)
   14f6b:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14f6f:	8b 40 48             	mov    0x48(%rax),%eax
   14f72:	89 45 a4             	mov    %eax,-0x5c(%rbp)
    for(int n=0; n<newSize; n++) {
   14f75:	8b 5d a0             	mov    -0x60(%rbp),%ebx
   14f78:	e8 c3 5a ff ff       	call   aa40 <omp_get_num_threads@plt>
   14f7d:	41 89 c4             	mov    %eax,%r12d
   14f80:	e8 fb 58 ff ff       	call   a880 <omp_get_thread_num@plt>
   14f85:	89 c6                	mov    %eax,%esi
   14f87:	89 d8                	mov    %ebx,%eax
   14f89:	99                   	cltd   
   14f8a:	41 f7 fc             	idiv   %r12d
   14f8d:	89 c1                	mov    %eax,%ecx
   14f8f:	89 d8                	mov    %ebx,%eax
   14f91:	99                   	cltd   
   14f92:	41 f7 fc             	idiv   %r12d
   14f95:	89 d0                	mov    %edx,%eax
   14f97:	39 c6                	cmp    %eax,%esi
   14f99:	0f 8c 07 02 00 00    	jl     151a6 <_Z7paddingPP10particle_tPii._omp_fn.2+0x2d7>
   14f9f:	89 f2                	mov    %esi,%edx
   14fa1:	0f af d1             	imul   %ecx,%edx
   14fa4:	01 d0                	add    %edx,%eax
   14fa6:	8d 14 08             	lea    (%rax,%rcx,1),%edx
   14fa9:	39 d0                	cmp    %edx,%eax
   14fab:	0f 8d 02 02 00 00    	jge    151b3 <_Z7paddingPP10particle_tPii._omp_fn.2+0x2e4>
   14fb1:	89 45 94             	mov    %eax,-0x6c(%rbp)
        if((newPos != oriPos) || (n>=oriSize))
   14fb4:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   14fb8:	48 3b 45 e0          	cmp    -0x20(%rbp),%rax
   14fbc:	0f 85 a3 01 00 00    	jne    15165 <_Z7paddingPP10particle_tPii._omp_fn.2+0x296>
   14fc2:	e9 d2 01 00 00       	jmp    15199 <_Z7paddingPP10particle_tPii._omp_fn.2+0x2ca>
        if((newVel != oriVel) || (n>=oriSize))
   14fc7:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14fcb:	48 3b 45 d8          	cmp    -0x28(%rbp),%rax
   14fcf:	0f 85 4f 01 00 00    	jne    15124 <_Z7paddingPP10particle_tPii._omp_fn.2+0x255>
   14fd5:	e9 7e 01 00 00       	jmp    15158 <_Z7paddingPP10particle_tPii._omp_fn.2+0x289>
        if((newAcc != oriAcc) || (n>=oriSize))
   14fda:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   14fde:	48 3b 45 d0          	cmp    -0x30(%rbp),%rax
   14fe2:	0f 85 fb 00 00 00    	jne    150e3 <_Z7paddingPP10particle_tPii._omp_fn.2+0x214>
   14fe8:	e9 2a 01 00 00       	jmp    15117 <_Z7paddingPP10particle_tPii._omp_fn.2+0x248>
        if((newFeat != oriFeat) || (n>=oriSize))
   14fed:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   14ff1:	48 3b 45 c8          	cmp    -0x38(%rbp),%rax
   14ff5:	0f 85 a7 00 00 00    	jne    150a2 <_Z7paddingPP10particle_tPii._omp_fn.2+0x1d3>
   14ffb:	e9 d6 00 00 00       	jmp    150d6 <_Z7paddingPP10particle_tPii._omp_fn.2+0x207>
        if(n>=oriSize) {
   15000:	8b 45 94             	mov    -0x6c(%rbp),%eax
   15003:	3b 45 a4             	cmp    -0x5c(%rbp),%eax
   15006:	7d 0e                	jge    15016 <_Z7paddingPP10particle_tPii._omp_fn.2+0x147>
   15008:	83 45 94 01          	addl   $0x1,-0x6c(%rbp)
   1500c:	39 55 94             	cmp    %edx,-0x6c(%rbp)
   1500f:	7c a3                	jl     14fb4 <_Z7paddingPP10particle_tPii._omp_fn.2+0xe5>
    #pragma omp parallel for
   15011:	e9 9d 01 00 00       	jmp    151b3 <_Z7paddingPP10particle_tPii._omp_fn.2+0x2e4>
            (newP+n)->enabled = false;
   15016:	8b 45 94             	mov    -0x6c(%rbp),%eax
   15019:	48 98                	cltq   
   1501b:	48 c1 e0 06          	shl    $0x6,%rax
   1501f:	48 89 c1             	mov    %rax,%rcx
   15022:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15026:	48 01 c8             	add    %rcx,%rax
   15029:	c6 40 38 00          	movb   $0x0,0x38(%rax)
            (newP+n)->id = UINT64_MAX - n;
   1502d:	8b 45 94             	mov    -0x6c(%rbp),%eax
   15030:	48 63 c8             	movslq %eax,%rcx
   15033:	8b 45 94             	mov    -0x6c(%rbp),%eax
   15036:	48 98                	cltq   
   15038:	48 c1 e0 06          	shl    $0x6,%rax
   1503c:	48 89 c6             	mov    %rax,%rsi
   1503f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15043:	48 01 f0             	add    %rsi,%rax
   15046:	48 f7 d1             	not    %rcx
   15049:	48 89 08             	mov    %rcx,(%rax)
            (newP+n)->ndim = ndim;
   1504c:	8b 45 94             	mov    -0x6c(%rbp),%eax
   1504f:	48 98                	cltq   
   15051:	48 c1 e0 06          	shl    $0x6,%rax
   15055:	48 89 c1             	mov    %rax,%rcx
   15058:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1505c:	48 01 c1             	add    %rax,%rcx
   1505f:	8b 45 9c             	mov    -0x64(%rbp),%eax
   15062:	89 41 30             	mov    %eax,0x30(%rcx)
            (newP+n)->nfeat = nfeat;
   15065:	8b 45 94             	mov    -0x6c(%rbp),%eax
   15068:	48 98                	cltq   
   1506a:	48 c1 e0 06          	shl    $0x6,%rax
   1506e:	48 89 c1             	mov    %rax,%rcx
   15071:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15075:	48 01 c1             	add    %rax,%rcx
   15078:	8b 45 98             	mov    -0x68(%rbp),%eax
   1507b:	89 41 34             	mov    %eax,0x34(%rcx)
            (newP+n)->updateAcceleration = newP->updateAcceleration;
   1507e:	8b 45 94             	mov    -0x6c(%rbp),%eax
   15081:	48 98                	cltq   
   15083:	48 c1 e0 06          	shl    $0x6,%rax
   15087:	48 89 c1             	mov    %rax,%rcx
   1508a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1508e:	48 01 c1             	add    %rax,%rcx
   15091:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15095:	48 8b 40 28          	mov    0x28(%rax),%rax
   15099:	48 89 41 28          	mov    %rax,0x28(%rcx)
   1509d:	e9 66 ff ff ff       	jmp    15008 <_Z7paddingPP10particle_tPii._omp_fn.2+0x139>
            (newP+n)->features = newFeat+n*nfeat;
   150a2:	8b 45 94             	mov    -0x6c(%rbp),%eax
   150a5:	0f af 45 98          	imul   -0x68(%rbp),%eax
   150a9:	48 98                	cltq   
   150ab:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   150b2:	00 
   150b3:	8b 45 94             	mov    -0x6c(%rbp),%eax
   150b6:	48 98                	cltq   
   150b8:	48 c1 e0 06          	shl    $0x6,%rax
   150bc:	48 89 c1             	mov    %rax,%rcx
   150bf:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   150c3:	48 01 c8             	add    %rcx,%rax
   150c6:	48 8b 4d a8          	mov    -0x58(%rbp),%rcx
   150ca:	48 01 f1             	add    %rsi,%rcx
   150cd:	48 89 48 20          	mov    %rcx,0x20(%rax)
   150d1:	e9 2a ff ff ff       	jmp    15000 <_Z7paddingPP10particle_tPii._omp_fn.2+0x131>
        if((newFeat != oriFeat) || (n>=oriSize))
   150d6:	8b 45 94             	mov    -0x6c(%rbp),%eax
   150d9:	3b 45 a4             	cmp    -0x5c(%rbp),%eax
   150dc:	7d c4                	jge    150a2 <_Z7paddingPP10particle_tPii._omp_fn.2+0x1d3>
   150de:	e9 1d ff ff ff       	jmp    15000 <_Z7paddingPP10particle_tPii._omp_fn.2+0x131>
            (newP+n)->acceleration = newAcc+n*ndim;
   150e3:	8b 45 94             	mov    -0x6c(%rbp),%eax
   150e6:	0f af 45 9c          	imul   -0x64(%rbp),%eax
   150ea:	48 98                	cltq   
   150ec:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   150f3:	00 
   150f4:	8b 45 94             	mov    -0x6c(%rbp),%eax
   150f7:	48 98                	cltq   
   150f9:	48 c1 e0 06          	shl    $0x6,%rax
   150fd:	48 89 c1             	mov    %rax,%rcx
   15100:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15104:	48 01 c8             	add    %rcx,%rax
   15107:	48 8b 4d b0          	mov    -0x50(%rbp),%rcx
   1510b:	48 01 f1             	add    %rsi,%rcx
   1510e:	48 89 48 18          	mov    %rcx,0x18(%rax)
   15112:	e9 d6 fe ff ff       	jmp    14fed <_Z7paddingPP10particle_tPii._omp_fn.2+0x11e>
        if((newAcc != oriAcc) || (n>=oriSize))
   15117:	8b 45 94             	mov    -0x6c(%rbp),%eax
   1511a:	3b 45 a4             	cmp    -0x5c(%rbp),%eax
   1511d:	7d c4                	jge    150e3 <_Z7paddingPP10particle_tPii._omp_fn.2+0x214>
   1511f:	e9 c9 fe ff ff       	jmp    14fed <_Z7paddingPP10particle_tPii._omp_fn.2+0x11e>
            (newP+n)->velocity = newVel+n*ndim;
   15124:	8b 45 94             	mov    -0x6c(%rbp),%eax
   15127:	0f af 45 9c          	imul   -0x64(%rbp),%eax
   1512b:	48 98                	cltq   
   1512d:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   15134:	00 
   15135:	8b 45 94             	mov    -0x6c(%rbp),%eax
   15138:	48 98                	cltq   
   1513a:	48 c1 e0 06          	shl    $0x6,%rax
   1513e:	48 89 c1             	mov    %rax,%rcx
   15141:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15145:	48 01 c8             	add    %rcx,%rax
   15148:	48 8b 4d b8          	mov    -0x48(%rbp),%rcx
   1514c:	48 01 f1             	add    %rsi,%rcx
   1514f:	48 89 48 10          	mov    %rcx,0x10(%rax)
   15153:	e9 82 fe ff ff       	jmp    14fda <_Z7paddingPP10particle_tPii._omp_fn.2+0x10b>
        if((newVel != oriVel) || (n>=oriSize))
   15158:	8b 45 94             	mov    -0x6c(%rbp),%eax
   1515b:	3b 45 a4             	cmp    -0x5c(%rbp),%eax
   1515e:	7d c4                	jge    15124 <_Z7paddingPP10particle_tPii._omp_fn.2+0x255>
   15160:	e9 75 fe ff ff       	jmp    14fda <_Z7paddingPP10particle_tPii._omp_fn.2+0x10b>
            (newP+n)->position = newPos+n*ndim;
   15165:	8b 45 94             	mov    -0x6c(%rbp),%eax
   15168:	0f af 45 9c          	imul   -0x64(%rbp),%eax
   1516c:	48 98                	cltq   
   1516e:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   15175:	00 
   15176:	8b 45 94             	mov    -0x6c(%rbp),%eax
   15179:	48 98                	cltq   
   1517b:	48 c1 e0 06          	shl    $0x6,%rax
   1517f:	48 89 c1             	mov    %rax,%rcx
   15182:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15186:	48 01 c8             	add    %rcx,%rax
   15189:	48 8b 4d c0          	mov    -0x40(%rbp),%rcx
   1518d:	48 01 f1             	add    %rsi,%rcx
   15190:	48 89 48 08          	mov    %rcx,0x8(%rax)
   15194:	e9 2e fe ff ff       	jmp    14fc7 <_Z7paddingPP10particle_tPii._omp_fn.2+0xf8>
        if((newPos != oriPos) || (n>=oriSize))
   15199:	8b 45 94             	mov    -0x6c(%rbp),%eax
   1519c:	3b 45 a4             	cmp    -0x5c(%rbp),%eax
   1519f:	7d c4                	jge    15165 <_Z7paddingPP10particle_tPii._omp_fn.2+0x296>
   151a1:	e9 21 fe ff ff       	jmp    14fc7 <_Z7paddingPP10particle_tPii._omp_fn.2+0xf8>
   151a6:	b8 00 00 00 00       	mov    $0x0,%eax
   151ab:	83 c1 01             	add    $0x1,%ecx
   151ae:	e9 ec fd ff ff       	jmp    14f9f <_Z7paddingPP10particle_tPii._omp_fn.2+0xd0>
    #pragma omp parallel for
   151b3:	90                   	nop
   151b4:	48 83 c4 70          	add    $0x70,%rsp
   151b8:	5b                   	pop    %rbx
   151b9:	41 5c                	pop    %r12
   151bb:	5d                   	pop    %rbp
   151bc:	c3                   	ret    

00000000000151bd <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.0>:
        #pragma omp parallel for
   151bd:	f3 0f 1e fa          	endbr64 
   151c1:	55                   	push   %rbp
   151c2:	48 89 e5             	mov    %rsp,%rbp
   151c5:	41 54                	push   %r12
   151c7:	53                   	push   %rbx
   151c8:	48 83 ec 50          	sub    $0x50,%rsp
   151cc:	48 89 7d a8          	mov    %rdi,-0x58(%rbp)
   151d0:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   151d4:	48 8b 40 28          	mov    0x28(%rax),%rax
   151d8:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   151dc:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   151e0:	48 8b 40 20          	mov    0x20(%rax),%rax
   151e4:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   151e8:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   151ec:	48 8b 40 18          	mov    0x18(%rax),%rax
   151f0:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   151f4:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   151f8:	48 8b 40 10          	mov    0x10(%rax),%rax
   151fc:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   15200:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15204:	8b 40 34             	mov    0x34(%rax),%eax
   15207:	89 45 b8             	mov    %eax,-0x48(%rbp)
   1520a:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   1520e:	8b 40 30             	mov    0x30(%rax),%eax
   15211:	89 45 bc             	mov    %eax,-0x44(%rbp)
   15214:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15218:	48 8b 40 08          	mov    0x8(%rax),%rax
   1521c:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   15220:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15224:	48 8b 00             	mov    (%rax),%rax
   15227:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
        for (int n=0; n<*localNParticlesX; n++) {
   1522b:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   1522f:	8b 18                	mov    (%rax),%ebx
   15231:	e8 0a 58 ff ff       	call   aa40 <omp_get_num_threads@plt>
   15236:	41 89 c4             	mov    %eax,%r12d
   15239:	e8 42 56 ff ff       	call   a880 <omp_get_thread_num@plt>
   1523e:	89 c6                	mov    %eax,%esi
   15240:	89 d8                	mov    %ebx,%eax
   15242:	99                   	cltd   
   15243:	41 f7 fc             	idiv   %r12d
   15246:	89 c1                	mov    %eax,%ecx
   15248:	89 d8                	mov    %ebx,%eax
   1524a:	99                   	cltd   
   1524b:	41 f7 fc             	idiv   %r12d
   1524e:	89 d0                	mov    %edx,%eax
   15250:	39 c6                	cmp    %eax,%esi
   15252:	0f 8c e0 00 00 00    	jl     15338 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.0+0x17b>
   15258:	89 f2                	mov    %esi,%edx
   1525a:	0f af d1             	imul   %ecx,%edx
   1525d:	01 d0                	add    %edx,%eax
   1525f:	8d 14 08             	lea    (%rax,%rcx,1),%edx
   15262:	39 d0                	cmp    %edx,%eax
   15264:	0f 8d db 00 00 00    	jge    15345 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.0+0x188>
   1526a:	89 45 b4             	mov    %eax,-0x4c(%rbp)
            ((*localParticlesX)+n)->position = pos+n*ndim;
   1526d:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15270:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   15274:	48 98                	cltq   
   15276:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   1527d:	00 
   1527e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15282:	48 8b 08             	mov    (%rax),%rcx
   15285:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15288:	48 98                	cltq   
   1528a:	48 c1 e0 06          	shl    $0x6,%rax
   1528e:	48 01 c8             	add    %rcx,%rax
   15291:	48 8b 4d d8          	mov    -0x28(%rbp),%rcx
   15295:	48 01 f1             	add    %rsi,%rcx
   15298:	48 89 48 08          	mov    %rcx,0x8(%rax)
            ((*localParticlesX)+n)->velocity = vel+n*ndim;
   1529c:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   1529f:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   152a3:	48 98                	cltq   
   152a5:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   152ac:	00 
   152ad:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   152b1:	48 8b 08             	mov    (%rax),%rcx
   152b4:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   152b7:	48 98                	cltq   
   152b9:	48 c1 e0 06          	shl    $0x6,%rax
   152bd:	48 01 c8             	add    %rcx,%rax
   152c0:	48 8b 4d d0          	mov    -0x30(%rbp),%rcx
   152c4:	48 01 f1             	add    %rsi,%rcx
   152c7:	48 89 48 10          	mov    %rcx,0x10(%rax)
            ((*localParticlesX)+n)->acceleration = acc+n*ndim;
   152cb:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   152ce:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   152d2:	48 98                	cltq   
   152d4:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   152db:	00 
   152dc:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   152e0:	48 8b 08             	mov    (%rax),%rcx
   152e3:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   152e6:	48 98                	cltq   
   152e8:	48 c1 e0 06          	shl    $0x6,%rax
   152ec:	48 01 c8             	add    %rcx,%rax
   152ef:	48 8b 4d c8          	mov    -0x38(%rbp),%rcx
   152f3:	48 01 f1             	add    %rsi,%rcx
   152f6:	48 89 48 18          	mov    %rcx,0x18(%rax)
            ((*localParticlesX)+n)->features = feat+n*nfeat;
   152fa:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   152fd:	0f af 45 b8          	imul   -0x48(%rbp),%eax
   15301:	48 98                	cltq   
   15303:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   1530a:	00 
   1530b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1530f:	48 8b 08             	mov    (%rax),%rcx
   15312:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15315:	48 98                	cltq   
   15317:	48 c1 e0 06          	shl    $0x6,%rax
   1531b:	48 01 c8             	add    %rcx,%rax
   1531e:	48 8b 4d c0          	mov    -0x40(%rbp),%rcx
   15322:	48 01 f1             	add    %rsi,%rcx
   15325:	48 89 48 20          	mov    %rcx,0x20(%rax)
   15329:	83 45 b4 01          	addl   $0x1,-0x4c(%rbp)
   1532d:	39 55 b4             	cmp    %edx,-0x4c(%rbp)
   15330:	0f 8c 37 ff ff ff    	jl     1526d <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.0+0xb0>
        #pragma omp parallel for
   15336:	eb 0d                	jmp    15345 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.0+0x188>
   15338:	b8 00 00 00 00       	mov    $0x0,%eax
   1533d:	83 c1 01             	add    $0x1,%ecx
   15340:	e9 13 ff ff ff       	jmp    15258 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.0+0x9b>
   15345:	90                   	nop
   15346:	48 83 c4 50          	add    $0x50,%rsp
   1534a:	5b                   	pop    %rbx
   1534b:	41 5c                	pop    %r12
   1534d:	5d                   	pop    %rbp
   1534e:	c3                   	ret    

000000000001534f <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.1>:
    #pragma omp parallel for
   1534f:	f3 0f 1e fa          	endbr64 
   15353:	55                   	push   %rbp
   15354:	48 89 e5             	mov    %rsp,%rbp
   15357:	41 54                	push   %r12
   15359:	53                   	push   %rbx
   1535a:	48 83 ec 50          	sub    $0x50,%rsp
   1535e:	48 89 7d a8          	mov    %rdi,-0x58(%rbp)
   15362:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15366:	48 8b 40 28          	mov    0x28(%rax),%rax
   1536a:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   1536e:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15372:	48 8b 40 20          	mov    0x20(%rax),%rax
   15376:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   1537a:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   1537e:	48 8b 40 18          	mov    0x18(%rax),%rax
   15382:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   15386:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   1538a:	48 8b 40 10          	mov    0x10(%rax),%rax
   1538e:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   15392:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15396:	8b 40 34             	mov    0x34(%rax),%eax
   15399:	89 45 b8             	mov    %eax,-0x48(%rbp)
   1539c:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   153a0:	8b 40 30             	mov    0x30(%rax),%eax
   153a3:	89 45 bc             	mov    %eax,-0x44(%rbp)
   153a6:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   153aa:	48 8b 40 08          	mov    0x8(%rax),%rax
   153ae:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   153b2:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   153b6:	48 8b 00             	mov    (%rax),%rax
   153b9:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    for (int n=0; n<*localNParticlesX; n++) {
   153bd:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   153c1:	8b 18                	mov    (%rax),%ebx
   153c3:	e8 78 56 ff ff       	call   aa40 <omp_get_num_threads@plt>
   153c8:	41 89 c4             	mov    %eax,%r12d
   153cb:	e8 b0 54 ff ff       	call   a880 <omp_get_thread_num@plt>
   153d0:	89 c6                	mov    %eax,%esi
   153d2:	89 d8                	mov    %ebx,%eax
   153d4:	99                   	cltd   
   153d5:	41 f7 fc             	idiv   %r12d
   153d8:	89 c1                	mov    %eax,%ecx
   153da:	89 d8                	mov    %ebx,%eax
   153dc:	99                   	cltd   
   153dd:	41 f7 fc             	idiv   %r12d
   153e0:	89 d0                	mov    %edx,%eax
   153e2:	39 c6                	cmp    %eax,%esi
   153e4:	0f 8c e0 00 00 00    	jl     154ca <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.1+0x17b>
   153ea:	89 f2                	mov    %esi,%edx
   153ec:	0f af d1             	imul   %ecx,%edx
   153ef:	01 d0                	add    %edx,%eax
   153f1:	8d 14 08             	lea    (%rax,%rcx,1),%edx
   153f4:	39 d0                	cmp    %edx,%eax
   153f6:	0f 8d db 00 00 00    	jge    154d7 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.1+0x188>
   153fc:	89 45 b4             	mov    %eax,-0x4c(%rbp)
        (*localParticlesX+n)->position = pos+n*ndim;
   153ff:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15402:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   15406:	48 98                	cltq   
   15408:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   1540f:	00 
   15410:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15414:	48 8b 08             	mov    (%rax),%rcx
   15417:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   1541a:	48 98                	cltq   
   1541c:	48 c1 e0 06          	shl    $0x6,%rax
   15420:	48 01 c8             	add    %rcx,%rax
   15423:	48 8b 4d d8          	mov    -0x28(%rbp),%rcx
   15427:	48 01 f1             	add    %rsi,%rcx
   1542a:	48 89 48 08          	mov    %rcx,0x8(%rax)
        (*localParticlesX+n)->velocity = vel+n*ndim;
   1542e:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15431:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   15435:	48 98                	cltq   
   15437:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   1543e:	00 
   1543f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15443:	48 8b 08             	mov    (%rax),%rcx
   15446:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15449:	48 98                	cltq   
   1544b:	48 c1 e0 06          	shl    $0x6,%rax
   1544f:	48 01 c8             	add    %rcx,%rax
   15452:	48 8b 4d d0          	mov    -0x30(%rbp),%rcx
   15456:	48 01 f1             	add    %rsi,%rcx
   15459:	48 89 48 10          	mov    %rcx,0x10(%rax)
        (*localParticlesX+n)->acceleration = acc+n*ndim;
   1545d:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15460:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   15464:	48 98                	cltq   
   15466:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   1546d:	00 
   1546e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15472:	48 8b 08             	mov    (%rax),%rcx
   15475:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15478:	48 98                	cltq   
   1547a:	48 c1 e0 06          	shl    $0x6,%rax
   1547e:	48 01 c8             	add    %rcx,%rax
   15481:	48 8b 4d c8          	mov    -0x38(%rbp),%rcx
   15485:	48 01 f1             	add    %rsi,%rcx
   15488:	48 89 48 18          	mov    %rcx,0x18(%rax)
        (*localParticlesX+n)->features = feat+n*nfeat;
   1548c:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   1548f:	0f af 45 b8          	imul   -0x48(%rbp),%eax
   15493:	48 98                	cltq   
   15495:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   1549c:	00 
   1549d:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   154a1:	48 8b 08             	mov    (%rax),%rcx
   154a4:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   154a7:	48 98                	cltq   
   154a9:	48 c1 e0 06          	shl    $0x6,%rax
   154ad:	48 01 c8             	add    %rcx,%rax
   154b0:	48 8b 4d c0          	mov    -0x40(%rbp),%rcx
   154b4:	48 01 f1             	add    %rsi,%rcx
   154b7:	48 89 48 20          	mov    %rcx,0x20(%rax)
   154bb:	83 45 b4 01          	addl   $0x1,-0x4c(%rbp)
   154bf:	39 55 b4             	cmp    %edx,-0x4c(%rbp)
   154c2:	0f 8c 37 ff ff ff    	jl     153ff <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.1+0xb0>
    #pragma omp parallel for
   154c8:	eb 0d                	jmp    154d7 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.1+0x188>
   154ca:	b8 00 00 00 00       	mov    $0x0,%eax
   154cf:	83 c1 01             	add    $0x1,%ecx
   154d2:	e9 13 ff ff ff       	jmp    153ea <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.1+0x9b>
   154d7:	90                   	nop
   154d8:	48 83 c4 50          	add    $0x50,%rsp
   154dc:	5b                   	pop    %rbx
   154dd:	41 5c                	pop    %r12
   154df:	5d                   	pop    %rbp
   154e0:	c3                   	ret    

00000000000154e1 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.2>:
        #pragma omp parallel for
   154e1:	f3 0f 1e fa          	endbr64 
   154e5:	55                   	push   %rbp
   154e6:	48 89 e5             	mov    %rsp,%rbp
   154e9:	41 54                	push   %r12
   154eb:	53                   	push   %rbx
   154ec:	48 83 ec 50          	sub    $0x50,%rsp
   154f0:	48 89 7d a8          	mov    %rdi,-0x58(%rbp)
   154f4:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   154f8:	48 8b 40 28          	mov    0x28(%rax),%rax
   154fc:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   15500:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15504:	48 8b 40 20          	mov    0x20(%rax),%rax
   15508:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   1550c:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15510:	48 8b 40 18          	mov    0x18(%rax),%rax
   15514:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   15518:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   1551c:	48 8b 40 10          	mov    0x10(%rax),%rax
   15520:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   15524:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15528:	8b 40 34             	mov    0x34(%rax),%eax
   1552b:	89 45 b8             	mov    %eax,-0x48(%rbp)
   1552e:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15532:	8b 40 30             	mov    0x30(%rax),%eax
   15535:	89 45 bc             	mov    %eax,-0x44(%rbp)
   15538:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   1553c:	48 8b 40 08          	mov    0x8(%rax),%rax
   15540:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   15544:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15548:	48 8b 00             	mov    (%rax),%rax
   1554b:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
        for (int n=0; n<*localNParticlesY; n++) {
   1554f:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   15553:	8b 18                	mov    (%rax),%ebx
   15555:	e8 e6 54 ff ff       	call   aa40 <omp_get_num_threads@plt>
   1555a:	41 89 c4             	mov    %eax,%r12d
   1555d:	e8 1e 53 ff ff       	call   a880 <omp_get_thread_num@plt>
   15562:	89 c6                	mov    %eax,%esi
   15564:	89 d8                	mov    %ebx,%eax
   15566:	99                   	cltd   
   15567:	41 f7 fc             	idiv   %r12d
   1556a:	89 c1                	mov    %eax,%ecx
   1556c:	89 d8                	mov    %ebx,%eax
   1556e:	99                   	cltd   
   1556f:	41 f7 fc             	idiv   %r12d
   15572:	89 d0                	mov    %edx,%eax
   15574:	39 c6                	cmp    %eax,%esi
   15576:	0f 8c e0 00 00 00    	jl     1565c <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.2+0x17b>
   1557c:	89 f2                	mov    %esi,%edx
   1557e:	0f af d1             	imul   %ecx,%edx
   15581:	01 d0                	add    %edx,%eax
   15583:	8d 14 08             	lea    (%rax,%rcx,1),%edx
   15586:	39 d0                	cmp    %edx,%eax
   15588:	0f 8d db 00 00 00    	jge    15669 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.2+0x188>
   1558e:	89 45 b4             	mov    %eax,-0x4c(%rbp)
            (*localParticlesY+n)->position = pos+n*ndim;
   15591:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15594:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   15598:	48 98                	cltq   
   1559a:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   155a1:	00 
   155a2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   155a6:	48 8b 08             	mov    (%rax),%rcx
   155a9:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   155ac:	48 98                	cltq   
   155ae:	48 c1 e0 06          	shl    $0x6,%rax
   155b2:	48 01 c8             	add    %rcx,%rax
   155b5:	48 8b 4d d8          	mov    -0x28(%rbp),%rcx
   155b9:	48 01 f1             	add    %rsi,%rcx
   155bc:	48 89 48 08          	mov    %rcx,0x8(%rax)
            (*localParticlesY+n)->velocity = vel+n*ndim;
   155c0:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   155c3:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   155c7:	48 98                	cltq   
   155c9:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   155d0:	00 
   155d1:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   155d5:	48 8b 08             	mov    (%rax),%rcx
   155d8:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   155db:	48 98                	cltq   
   155dd:	48 c1 e0 06          	shl    $0x6,%rax
   155e1:	48 01 c8             	add    %rcx,%rax
   155e4:	48 8b 4d d0          	mov    -0x30(%rbp),%rcx
   155e8:	48 01 f1             	add    %rsi,%rcx
   155eb:	48 89 48 10          	mov    %rcx,0x10(%rax)
            (*localParticlesY+n)->acceleration = acc+n*ndim;
   155ef:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   155f2:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   155f6:	48 98                	cltq   
   155f8:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   155ff:	00 
   15600:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15604:	48 8b 08             	mov    (%rax),%rcx
   15607:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   1560a:	48 98                	cltq   
   1560c:	48 c1 e0 06          	shl    $0x6,%rax
   15610:	48 01 c8             	add    %rcx,%rax
   15613:	48 8b 4d c8          	mov    -0x38(%rbp),%rcx
   15617:	48 01 f1             	add    %rsi,%rcx
   1561a:	48 89 48 18          	mov    %rcx,0x18(%rax)
            (*localParticlesY+n)->features = feat+n*nfeat;
   1561e:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15621:	0f af 45 b8          	imul   -0x48(%rbp),%eax
   15625:	48 98                	cltq   
   15627:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   1562e:	00 
   1562f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15633:	48 8b 08             	mov    (%rax),%rcx
   15636:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15639:	48 98                	cltq   
   1563b:	48 c1 e0 06          	shl    $0x6,%rax
   1563f:	48 01 c8             	add    %rcx,%rax
   15642:	48 8b 4d c0          	mov    -0x40(%rbp),%rcx
   15646:	48 01 f1             	add    %rsi,%rcx
   15649:	48 89 48 20          	mov    %rcx,0x20(%rax)
   1564d:	83 45 b4 01          	addl   $0x1,-0x4c(%rbp)
   15651:	39 55 b4             	cmp    %edx,-0x4c(%rbp)
   15654:	0f 8c 37 ff ff ff    	jl     15591 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.2+0xb0>
        #pragma omp parallel for
   1565a:	eb 0d                	jmp    15669 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.2+0x188>
   1565c:	b8 00 00 00 00       	mov    $0x0,%eax
   15661:	83 c1 01             	add    $0x1,%ecx
   15664:	e9 13 ff ff ff       	jmp    1557c <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.2+0x9b>
   15669:	90                   	nop
   1566a:	48 83 c4 50          	add    $0x50,%rsp
   1566e:	5b                   	pop    %rbx
   1566f:	41 5c                	pop    %r12
   15671:	5d                   	pop    %rbp
   15672:	c3                   	ret    

0000000000015673 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.3>:
    #pragma omp parallel for
   15673:	f3 0f 1e fa          	endbr64 
   15677:	55                   	push   %rbp
   15678:	48 89 e5             	mov    %rsp,%rbp
   1567b:	41 54                	push   %r12
   1567d:	53                   	push   %rbx
   1567e:	48 83 ec 50          	sub    $0x50,%rsp
   15682:	48 89 7d a8          	mov    %rdi,-0x58(%rbp)
   15686:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   1568a:	48 8b 40 28          	mov    0x28(%rax),%rax
   1568e:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   15692:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15696:	48 8b 40 20          	mov    0x20(%rax),%rax
   1569a:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   1569e:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   156a2:	48 8b 40 18          	mov    0x18(%rax),%rax
   156a6:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   156aa:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   156ae:	48 8b 40 10          	mov    0x10(%rax),%rax
   156b2:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   156b6:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   156ba:	8b 40 34             	mov    0x34(%rax),%eax
   156bd:	89 45 b8             	mov    %eax,-0x48(%rbp)
   156c0:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   156c4:	8b 40 30             	mov    0x30(%rax),%eax
   156c7:	89 45 bc             	mov    %eax,-0x44(%rbp)
   156ca:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   156ce:	48 8b 40 08          	mov    0x8(%rax),%rax
   156d2:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   156d6:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   156da:	48 8b 00             	mov    (%rax),%rax
   156dd:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    for (int n=0; n<*localNParticlesY; n++) {
   156e1:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   156e5:	8b 18                	mov    (%rax),%ebx
   156e7:	e8 54 53 ff ff       	call   aa40 <omp_get_num_threads@plt>
   156ec:	41 89 c4             	mov    %eax,%r12d
   156ef:	e8 8c 51 ff ff       	call   a880 <omp_get_thread_num@plt>
   156f4:	89 c6                	mov    %eax,%esi
   156f6:	89 d8                	mov    %ebx,%eax
   156f8:	99                   	cltd   
   156f9:	41 f7 fc             	idiv   %r12d
   156fc:	89 c1                	mov    %eax,%ecx
   156fe:	89 d8                	mov    %ebx,%eax
   15700:	99                   	cltd   
   15701:	41 f7 fc             	idiv   %r12d
   15704:	89 d0                	mov    %edx,%eax
   15706:	39 c6                	cmp    %eax,%esi
   15708:	0f 8c e0 00 00 00    	jl     157ee <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.3+0x17b>
   1570e:	89 f2                	mov    %esi,%edx
   15710:	0f af d1             	imul   %ecx,%edx
   15713:	01 d0                	add    %edx,%eax
   15715:	8d 14 08             	lea    (%rax,%rcx,1),%edx
   15718:	39 d0                	cmp    %edx,%eax
   1571a:	0f 8d db 00 00 00    	jge    157fb <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.3+0x188>
   15720:	89 45 b4             	mov    %eax,-0x4c(%rbp)
        (*localParticlesY+n)->position = pos+n*ndim;
   15723:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15726:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   1572a:	48 98                	cltq   
   1572c:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   15733:	00 
   15734:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15738:	48 8b 08             	mov    (%rax),%rcx
   1573b:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   1573e:	48 98                	cltq   
   15740:	48 c1 e0 06          	shl    $0x6,%rax
   15744:	48 01 c8             	add    %rcx,%rax
   15747:	48 8b 4d d8          	mov    -0x28(%rbp),%rcx
   1574b:	48 01 f1             	add    %rsi,%rcx
   1574e:	48 89 48 08          	mov    %rcx,0x8(%rax)
        (*localParticlesY+n)->velocity = vel+n*ndim;
   15752:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15755:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   15759:	48 98                	cltq   
   1575b:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   15762:	00 
   15763:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15767:	48 8b 08             	mov    (%rax),%rcx
   1576a:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   1576d:	48 98                	cltq   
   1576f:	48 c1 e0 06          	shl    $0x6,%rax
   15773:	48 01 c8             	add    %rcx,%rax
   15776:	48 8b 4d d0          	mov    -0x30(%rbp),%rcx
   1577a:	48 01 f1             	add    %rsi,%rcx
   1577d:	48 89 48 10          	mov    %rcx,0x10(%rax)
        (*localParticlesY+n)->acceleration = acc+n*ndim;
   15781:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   15784:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   15788:	48 98                	cltq   
   1578a:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   15791:	00 
   15792:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15796:	48 8b 08             	mov    (%rax),%rcx
   15799:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   1579c:	48 98                	cltq   
   1579e:	48 c1 e0 06          	shl    $0x6,%rax
   157a2:	48 01 c8             	add    %rcx,%rax
   157a5:	48 8b 4d c8          	mov    -0x38(%rbp),%rcx
   157a9:	48 01 f1             	add    %rsi,%rcx
   157ac:	48 89 48 18          	mov    %rcx,0x18(%rax)
        (*localParticlesY+n)->features = feat+n*nfeat;
   157b0:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   157b3:	0f af 45 b8          	imul   -0x48(%rbp),%eax
   157b7:	48 98                	cltq   
   157b9:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   157c0:	00 
   157c1:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   157c5:	48 8b 08             	mov    (%rax),%rcx
   157c8:	8b 45 b4             	mov    -0x4c(%rbp),%eax
   157cb:	48 98                	cltq   
   157cd:	48 c1 e0 06          	shl    $0x6,%rax
   157d1:	48 01 c8             	add    %rcx,%rax
   157d4:	48 8b 4d c0          	mov    -0x40(%rbp),%rcx
   157d8:	48 01 f1             	add    %rsi,%rcx
   157db:	48 89 48 20          	mov    %rcx,0x20(%rax)
   157df:	83 45 b4 01          	addl   $0x1,-0x4c(%rbp)
   157e3:	39 55 b4             	cmp    %edx,-0x4c(%rbp)
   157e6:	0f 8c 37 ff ff ff    	jl     15723 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.3+0xb0>
    #pragma omp parallel for
   157ec:	eb 0d                	jmp    157fb <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.3+0x188>
   157ee:	b8 00 00 00 00       	mov    $0x0,%eax
   157f3:	83 c1 01             	add    $0x1,%ecx
   157f6:	e9 13 ff ff ff       	jmp    1570e <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t._omp_fn.3+0x9b>
   157fb:	90                   	nop
   157fc:	48 83 c4 50          	add    $0x50,%rsp
   15800:	5b                   	pop    %rbx
   15801:	41 5c                	pop    %r12
   15803:	5d                   	pop    %rbp
   15804:	c3                   	ret    

0000000000015805 <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.0>:
    #pragma omp parallel for collapse(2)
   15805:	f3 0f 1e fa          	endbr64 
   15809:	55                   	push   %rbp
   1580a:	48 89 e5             	mov    %rsp,%rbp
   1580d:	41 56                	push   %r14
   1580f:	41 55                	push   %r13
   15811:	41 54                	push   %r12
   15813:	53                   	push   %rbx
   15814:	48 83 ec 30          	sub    $0x30,%rsp
   15818:	48 89 7d b8          	mov    %rdi,-0x48(%rbp)
   1581c:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   15820:	48 8b 40 08          	mov    0x8(%rax),%rax
   15824:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   15828:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   1582c:	8b 40 10             	mov    0x10(%rax),%eax
   1582f:	89 45 cc             	mov    %eax,-0x34(%rbp)
   15832:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   15836:	48 8b 00             	mov    (%rax),%rax
   15839:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
        for(int nn=0; nn<localParticlesX->ndim; nn++)
   1583d:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   15841:	8b 58 30             	mov    0x30(%rax),%ebx
    for(int n=0; n<localNParticlesX; n++)
   15844:	8b 45 cc             	mov    -0x34(%rbp),%eax
   15847:	85 c0                	test   %eax,%eax
   15849:	0f 8e ed 00 00 00    	jle    1593c <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.0+0x137>
   1584f:	85 db                	test   %ebx,%ebx
   15851:	0f 8e e5 00 00 00    	jle    1593c <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.0+0x137>
   15857:	41 89 dd             	mov    %ebx,%r13d
   1585a:	41 0f af c5          	imul   %r13d,%eax
   1585e:	41 89 c4             	mov    %eax,%r12d
   15861:	e8 da 51 ff ff       	call   aa40 <omp_get_num_threads@plt>
   15866:	41 89 c6             	mov    %eax,%r14d
   15869:	e8 12 50 ff ff       	call   a880 <omp_get_thread_num@plt>
   1586e:	89 c1                	mov    %eax,%ecx
   15870:	44 89 e0             	mov    %r12d,%eax
   15873:	ba 00 00 00 00       	mov    $0x0,%edx
   15878:	41 f7 f6             	div    %r14d
   1587b:	89 c6                	mov    %eax,%esi
   1587d:	44 89 e0             	mov    %r12d,%eax
   15880:	ba 00 00 00 00       	mov    $0x0,%edx
   15885:	41 f7 f6             	div    %r14d
   15888:	89 d0                	mov    %edx,%eax
   1588a:	39 c1                	cmp    %eax,%ecx
   1588c:	0f 82 9d 00 00 00    	jb     1592f <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.0+0x12a>
   15892:	0f af ce             	imul   %esi,%ecx
   15895:	89 ca                	mov    %ecx,%edx
   15897:	8d 0c 02             	lea    (%rdx,%rax,1),%ecx
   1589a:	01 ce                	add    %ecx,%esi
   1589c:	39 f1                	cmp    %esi,%ecx
   1589e:	0f 83 98 00 00 00    	jae    1593c <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.0+0x137>
   158a4:	89 c8                	mov    %ecx,%eax
   158a6:	ba 00 00 00 00       	mov    $0x0,%edx
   158ab:	41 f7 f5             	div    %r13d
   158ae:	89 d0                	mov    %edx,%eax
   158b0:	89 45 c4             	mov    %eax,-0x3c(%rbp)
   158b3:	89 c8                	mov    %ecx,%eax
   158b5:	ba 00 00 00 00       	mov    $0x0,%edx
   158ba:	41 f7 f5             	div    %r13d
   158bd:	89 45 c8             	mov    %eax,-0x38(%rbp)
                    *(recvBuf1 + n*localParticlesX->ndim + nn);
   158c0:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   158c4:	8b 40 30             	mov    0x30(%rax),%eax
   158c7:	0f af 45 c8          	imul   -0x38(%rbp),%eax
   158cb:	48 63 d0             	movslq %eax,%rdx
   158ce:	8b 45 c4             	mov    -0x3c(%rbp),%eax
   158d1:	48 98                	cltq   
   158d3:	48 01 d0             	add    %rdx,%rax
   158d6:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   158dd:	00 
   158de:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   158e2:	48 01 c2             	add    %rax,%rdx
            *((localParticlesX->acceleration) + n*localParticlesX->ndim + nn) = 
   158e5:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   158e9:	48 8b 78 18          	mov    0x18(%rax),%rdi
   158ed:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   158f1:	8b 40 30             	mov    0x30(%rax),%eax
   158f4:	0f af 45 c8          	imul   -0x38(%rbp),%eax
   158f8:	4c 63 c0             	movslq %eax,%r8
   158fb:	8b 45 c4             	mov    -0x3c(%rbp),%eax
   158fe:	48 98                	cltq   
   15900:	4c 01 c0             	add    %r8,%rax
   15903:	48 c1 e0 03          	shl    $0x3,%rax
   15907:	48 01 f8             	add    %rdi,%rax
                    *(recvBuf1 + n*localParticlesX->ndim + nn);
   1590a:	f2 0f 10 02          	movsd  (%rdx),%xmm0
            *((localParticlesX->acceleration) + n*localParticlesX->ndim + nn) = 
   1590e:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   15912:	83 c1 01             	add    $0x1,%ecx
   15915:	39 f1                	cmp    %esi,%ecx
   15917:	73 23                	jae    1593c <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.0+0x137>
   15919:	83 45 c4 01          	addl   $0x1,-0x3c(%rbp)
   1591d:	39 5d c4             	cmp    %ebx,-0x3c(%rbp)
   15920:	7c 9e                	jl     158c0 <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.0+0xbb>
   15922:	c7 45 c4 00 00 00 00 	movl   $0x0,-0x3c(%rbp)
   15929:	83 45 c8 01          	addl   $0x1,-0x38(%rbp)
   1592d:	eb 91                	jmp    158c0 <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.0+0xbb>
   1592f:	b8 00 00 00 00       	mov    $0x0,%eax
   15934:	83 c6 01             	add    $0x1,%esi
   15937:	e9 56 ff ff ff       	jmp    15892 <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.0+0x8d>
    #pragma omp parallel for collapse(2)
   1593c:	90                   	nop
   1593d:	48 83 c4 30          	add    $0x30,%rsp
   15941:	5b                   	pop    %rbx
   15942:	41 5c                	pop    %r12
   15944:	41 5d                	pop    %r13
   15946:	41 5e                	pop    %r14
   15948:	5d                   	pop    %rbp
   15949:	c3                   	ret    

000000000001594a <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.1>:
    #pragma omp parallel for collapse(2)
   1594a:	f3 0f 1e fa          	endbr64 
   1594e:	55                   	push   %rbp
   1594f:	48 89 e5             	mov    %rsp,%rbp
   15952:	41 56                	push   %r14
   15954:	41 55                	push   %r13
   15956:	41 54                	push   %r12
   15958:	53                   	push   %rbx
   15959:	48 83 ec 30          	sub    $0x30,%rsp
   1595d:	48 89 7d b8          	mov    %rdi,-0x48(%rbp)
   15961:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   15965:	48 8b 40 08          	mov    0x8(%rax),%rax
   15969:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   1596d:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   15971:	8b 40 10             	mov    0x10(%rax),%eax
   15974:	89 45 cc             	mov    %eax,-0x34(%rbp)
   15977:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   1597b:	48 8b 00             	mov    (%rax),%rax
   1597e:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
        for(int nn=0; nn<localParticlesY->ndim; nn++)
   15982:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   15986:	8b 58 30             	mov    0x30(%rax),%ebx
    for(int n=0; n<localNParticlesY; n++)
   15989:	8b 45 cc             	mov    -0x34(%rbp),%eax
   1598c:	85 c0                	test   %eax,%eax
   1598e:	0f 8e ed 00 00 00    	jle    15a81 <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.1+0x137>
   15994:	85 db                	test   %ebx,%ebx
   15996:	0f 8e e5 00 00 00    	jle    15a81 <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.1+0x137>
   1599c:	41 89 dd             	mov    %ebx,%r13d
   1599f:	41 0f af c5          	imul   %r13d,%eax
   159a3:	41 89 c4             	mov    %eax,%r12d
   159a6:	e8 95 50 ff ff       	call   aa40 <omp_get_num_threads@plt>
   159ab:	41 89 c6             	mov    %eax,%r14d
   159ae:	e8 cd 4e ff ff       	call   a880 <omp_get_thread_num@plt>
   159b3:	89 c1                	mov    %eax,%ecx
   159b5:	44 89 e0             	mov    %r12d,%eax
   159b8:	ba 00 00 00 00       	mov    $0x0,%edx
   159bd:	41 f7 f6             	div    %r14d
   159c0:	89 c6                	mov    %eax,%esi
   159c2:	44 89 e0             	mov    %r12d,%eax
   159c5:	ba 00 00 00 00       	mov    $0x0,%edx
   159ca:	41 f7 f6             	div    %r14d
   159cd:	89 d0                	mov    %edx,%eax
   159cf:	39 c1                	cmp    %eax,%ecx
   159d1:	0f 82 9d 00 00 00    	jb     15a74 <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.1+0x12a>
   159d7:	0f af ce             	imul   %esi,%ecx
   159da:	89 ca                	mov    %ecx,%edx
   159dc:	8d 0c 02             	lea    (%rdx,%rax,1),%ecx
   159df:	01 ce                	add    %ecx,%esi
   159e1:	39 f1                	cmp    %esi,%ecx
   159e3:	0f 83 98 00 00 00    	jae    15a81 <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.1+0x137>
   159e9:	89 c8                	mov    %ecx,%eax
   159eb:	ba 00 00 00 00       	mov    $0x0,%edx
   159f0:	41 f7 f5             	div    %r13d
   159f3:	89 d0                	mov    %edx,%eax
   159f5:	89 45 c4             	mov    %eax,-0x3c(%rbp)
   159f8:	89 c8                	mov    %ecx,%eax
   159fa:	ba 00 00 00 00       	mov    $0x0,%edx
   159ff:	41 f7 f5             	div    %r13d
   15a02:	89 45 c8             	mov    %eax,-0x38(%rbp)
                    *(recvBuf2+n*localParticlesY->ndim+nn);
   15a05:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   15a09:	8b 40 30             	mov    0x30(%rax),%eax
   15a0c:	0f af 45 c8          	imul   -0x38(%rbp),%eax
   15a10:	48 63 d0             	movslq %eax,%rdx
   15a13:	8b 45 c4             	mov    -0x3c(%rbp),%eax
   15a16:	48 98                	cltq   
   15a18:	48 01 d0             	add    %rdx,%rax
   15a1b:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   15a22:	00 
   15a23:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   15a27:	48 01 c2             	add    %rax,%rdx
            *((localParticlesY->acceleration)+n*localParticlesY->ndim+nn) = 
   15a2a:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   15a2e:	48 8b 78 18          	mov    0x18(%rax),%rdi
   15a32:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   15a36:	8b 40 30             	mov    0x30(%rax),%eax
   15a39:	0f af 45 c8          	imul   -0x38(%rbp),%eax
   15a3d:	4c 63 c0             	movslq %eax,%r8
   15a40:	8b 45 c4             	mov    -0x3c(%rbp),%eax
   15a43:	48 98                	cltq   
   15a45:	4c 01 c0             	add    %r8,%rax
   15a48:	48 c1 e0 03          	shl    $0x3,%rax
   15a4c:	48 01 f8             	add    %rdi,%rax
                    *(recvBuf2+n*localParticlesY->ndim+nn);
   15a4f:	f2 0f 10 02          	movsd  (%rdx),%xmm0
            *((localParticlesY->acceleration)+n*localParticlesY->ndim+nn) = 
   15a53:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   15a57:	83 c1 01             	add    $0x1,%ecx
   15a5a:	39 f1                	cmp    %esi,%ecx
   15a5c:	73 23                	jae    15a81 <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.1+0x137>
   15a5e:	83 45 c4 01          	addl   $0x1,-0x3c(%rbp)
   15a62:	39 5d c4             	cmp    %ebx,-0x3c(%rbp)
   15a65:	7c 9e                	jl     15a05 <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.1+0xbb>
   15a67:	c7 45 c4 00 00 00 00 	movl   $0x0,-0x3c(%rbp)
   15a6e:	83 45 c8 01          	addl   $0x1,-0x38(%rbp)
   15a72:	eb 91                	jmp    15a05 <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.1+0xbb>
   15a74:	b8 00 00 00 00       	mov    $0x0,%eax
   15a79:	83 c6 01             	add    $0x1,%esi
   15a7c:	e9 56 ff ff ff       	jmp    159d7 <_Z6reduceP10particle_tiS0_iP10topology_t._omp_fn.1+0x8d>
    #pragma omp parallel for collapse(2)
   15a81:	90                   	nop
   15a82:	48 83 c4 30          	add    $0x30,%rsp
   15a86:	5b                   	pop    %rbx
   15a87:	41 5c                	pop    %r12
   15a89:	41 5d                	pop    %r13
   15a8b:	41 5e                	pop    %r14
   15a8d:	5d                   	pop    %rbp
   15a8e:	c3                   	ret    

0000000000015a8f <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t._omp_fn.0>:
        #pragma omp parallel for
   15a8f:	f3 0f 1e fa          	endbr64 
   15a93:	55                   	push   %rbp
   15a94:	48 89 e5             	mov    %rsp,%rbp
   15a97:	41 54                	push   %r12
   15a99:	53                   	push   %rbx
   15a9a:	48 83 ec 50          	sub    $0x50,%rsp
   15a9e:	48 89 7d a8          	mov    %rdi,-0x58(%rbp)
   15aa2:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15aa6:	48 8b 40 28          	mov    0x28(%rax),%rax
   15aaa:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   15aae:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15ab2:	48 8b 40 20          	mov    0x20(%rax),%rax
   15ab6:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   15aba:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15abe:	48 8b 40 18          	mov    0x18(%rax),%rax
   15ac2:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   15ac6:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15aca:	48 8b 40 10          	mov    0x10(%rax),%rax
   15ace:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   15ad2:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15ad6:	48 8b 40 08          	mov    0x8(%rax),%rax
   15ada:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   15ade:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   15ae2:	48 8b 00             	mov    (%rax),%rax
   15ae5:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
        for(int n=0; n<(*globalNParticles); n++) {
   15ae9:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   15aed:	8b 18                	mov    (%rax),%ebx
   15aef:	e8 4c 4f ff ff       	call   aa40 <omp_get_num_threads@plt>
   15af4:	41 89 c4             	mov    %eax,%r12d
   15af7:	e8 84 4d ff ff       	call   a880 <omp_get_thread_num@plt>
   15afc:	89 c6                	mov    %eax,%esi
   15afe:	89 d8                	mov    %ebx,%eax
   15b00:	99                   	cltd   
   15b01:	41 f7 fc             	idiv   %r12d
   15b04:	89 c1                	mov    %eax,%ecx
   15b06:	89 d8                	mov    %ebx,%eax
   15b08:	99                   	cltd   
   15b09:	41 f7 fc             	idiv   %r12d
   15b0c:	89 d0                	mov    %edx,%eax
   15b0e:	39 c6                	cmp    %eax,%esi
   15b10:	0f 8c cc 00 00 00    	jl     15be2 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t._omp_fn.0+0x153>
   15b16:	89 f2                	mov    %esi,%edx
   15b18:	0f af d1             	imul   %ecx,%edx
   15b1b:	01 d0                	add    %edx,%eax
   15b1d:	8d 14 08             	lea    (%rax,%rcx,1),%edx
   15b20:	39 d0                	cmp    %edx,%eax
   15b22:	0f 8d c7 00 00 00    	jge    15bef <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t._omp_fn.0+0x160>
   15b28:	89 45 bc             	mov    %eax,-0x44(%rbp)
            (*globalParticles)->position = pos+n*((*globalParticles)->ndim);
   15b2b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15b2f:	48 8b 00             	mov    (%rax),%rax
   15b32:	8b 40 30             	mov    0x30(%rax),%eax
   15b35:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   15b39:	48 98                	cltq   
   15b3b:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   15b42:	00 
   15b43:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15b47:	48 8b 00             	mov    (%rax),%rax
   15b4a:	48 8b 4d d8          	mov    -0x28(%rbp),%rcx
   15b4e:	48 01 f1             	add    %rsi,%rcx
   15b51:	48 89 48 08          	mov    %rcx,0x8(%rax)
            (*globalParticles)->velocity = vel+n*((*globalParticles)->ndim);
   15b55:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15b59:	48 8b 00             	mov    (%rax),%rax
   15b5c:	8b 40 30             	mov    0x30(%rax),%eax
   15b5f:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   15b63:	48 98                	cltq   
   15b65:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   15b6c:	00 
   15b6d:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15b71:	48 8b 00             	mov    (%rax),%rax
   15b74:	48 8b 4d d0          	mov    -0x30(%rbp),%rcx
   15b78:	48 01 f1             	add    %rsi,%rcx
   15b7b:	48 89 48 10          	mov    %rcx,0x10(%rax)
            (*globalParticles)->acceleration = acc+n*((*globalParticles)->ndim);
   15b7f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15b83:	48 8b 00             	mov    (%rax),%rax
   15b86:	8b 40 30             	mov    0x30(%rax),%eax
   15b89:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   15b8d:	48 98                	cltq   
   15b8f:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   15b96:	00 
   15b97:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15b9b:	48 8b 00             	mov    (%rax),%rax
   15b9e:	48 8b 4d c8          	mov    -0x38(%rbp),%rcx
   15ba2:	48 01 f1             	add    %rsi,%rcx
   15ba5:	48 89 48 18          	mov    %rcx,0x18(%rax)
            (*globalParticles)->features = feat+n*((*globalParticles)->nfeat);
   15ba9:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15bad:	48 8b 00             	mov    (%rax),%rax
   15bb0:	8b 40 34             	mov    0x34(%rax),%eax
   15bb3:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   15bb7:	48 98                	cltq   
   15bb9:	48 8d 34 c5 00 00 00 	lea    0x0(,%rax,8),%rsi
   15bc0:	00 
   15bc1:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   15bc5:	48 8b 00             	mov    (%rax),%rax
   15bc8:	48 8b 4d c0          	mov    -0x40(%rbp),%rcx
   15bcc:	48 01 f1             	add    %rsi,%rcx
   15bcf:	48 89 48 20          	mov    %rcx,0x20(%rax)
   15bd3:	83 45 bc 01          	addl   $0x1,-0x44(%rbp)
   15bd7:	39 55 bc             	cmp    %edx,-0x44(%rbp)
   15bda:	0f 8c 4b ff ff ff    	jl     15b2b <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t._omp_fn.0+0x9c>
        #pragma omp parallel for
   15be0:	eb 0d                	jmp    15bef <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t._omp_fn.0+0x160>
   15be2:	b8 00 00 00 00       	mov    $0x0,%eax
   15be7:	83 c1 01             	add    $0x1,%ecx
   15bea:	e9 27 ff ff ff       	jmp    15b16 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t._omp_fn.0+0x87>
   15bef:	90                   	nop
   15bf0:	48 83 c4 50          	add    $0x50,%rsp
   15bf4:	5b                   	pop    %rbx
   15bf5:	41 5c                	pop    %r12
   15bf7:	5d                   	pop    %rbp
   15bf8:	c3                   	ret    

0000000000015bf9 <_Z16sizeAfterPaddingii>:
__inline__ int sizeAfterPadding(int nParticles, int minimalGrid) {
   15bf9:	f3 0f 1e fa          	endbr64 
   15bfd:	55                   	push   %rbp
   15bfe:	48 89 e5             	mov    %rsp,%rbp
   15c01:	89 7d fc             	mov    %edi,-0x4(%rbp)
   15c04:	89 75 f8             	mov    %esi,-0x8(%rbp)
    return (nParticles+minimalGrid-1)/minimalGrid * minimalGrid;
   15c07:	8b 55 fc             	mov    -0x4(%rbp),%edx
   15c0a:	8b 45 f8             	mov    -0x8(%rbp),%eax
   15c0d:	01 d0                	add    %edx,%eax
   15c0f:	83 e8 01             	sub    $0x1,%eax
   15c12:	99                   	cltd   
   15c13:	f7 7d f8             	idivl  -0x8(%rbp)
   15c16:	0f af 45 f8          	imul   -0x8(%rbp),%eax
}
   15c1a:	5d                   	pop    %rbp
   15c1b:	c3                   	ret    

Disassembly of section .fini:

0000000000015c1c <_fini>:
   15c1c:	f3 0f 1e fa          	endbr64 
   15c20:	48 83 ec 08          	sub    $0x8,%rsp
   15c24:	48 83 c4 08          	add    $0x8,%rsp
   15c28:	c3                   	ret    
