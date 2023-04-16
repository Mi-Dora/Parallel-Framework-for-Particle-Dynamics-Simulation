
build/test/TEST_dense_mpi:     file format elf64-x86-64


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
    9020:	ff 35 da 38 01 00    	push   0x138da(%rip)        # 1c900 <_GLOBAL_OFFSET_TABLE_+0x8>
    9026:	f2 ff 25 db 38 01 00 	bnd jmp *0x138db(%rip)        # 1c908 <_GLOBAL_OFFSET_TABLE_+0x10>
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

Disassembly of section .plt.got:

0000000000009d90 <__cxa_finalize@plt>:
    9d90:	f3 0f 1e fa          	endbr64 
    9d94:	f2 ff 25 2d 32 01 00 	bnd jmp *0x1322d(%rip)        # 1cfc8 <__cxa_finalize@GLIBC_2.2.5>
    9d9b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

Disassembly of section .plt.sec:

0000000000009da0 <_Znam@plt>:
    9da0:	f3 0f 1e fa          	endbr64 
    9da4:	f2 ff 25 65 2b 01 00 	bnd jmp *0x12b65(%rip)        # 1c910 <_Znam@GLIBCXX_3.4>
    9dab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009db0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE8capacityEv@plt>:
    9db0:	f3 0f 1e fa          	endbr64 
    9db4:	f2 ff 25 5d 2b 01 00 	bnd jmp *0x12b5d(%rip)        # 1c918 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE8capacityEv@GLIBCXX_3.4.21>
    9dbb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009dc0 <printf@plt>:
    9dc0:	f3 0f 1e fa          	endbr64 
    9dc4:	f2 ff 25 55 2b 01 00 	bnd jmp *0x12b55(%rip)        # 1c920 <printf@GLIBC_2.2.5>
    9dcb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009dd0 <MPI_Win_unlock@plt>:
    9dd0:	f3 0f 1e fa          	endbr64 
    9dd4:	f2 ff 25 4d 2b 01 00 	bnd jmp *0x12b4d(%rip)        # 1c928 <MPI_Win_unlock@Base>
    9ddb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009de0 <MPI_Group_excl@plt>:
    9de0:	f3 0f 1e fa          	endbr64 
    9de4:	f2 ff 25 45 2b 01 00 	bnd jmp *0x12b45(%rip)        # 1c930 <MPI_Group_excl@Base>
    9deb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009df0 <MPI_Allgather@plt>:
    9df0:	f3 0f 1e fa          	endbr64 
    9df4:	f2 ff 25 3d 2b 01 00 	bnd jmp *0x12b3d(%rip)        # 1c938 <MPI_Allgather@Base>
    9dfb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e00 <MPI_Recv_init@plt>:
    9e00:	f3 0f 1e fa          	endbr64 
    9e04:	f2 ff 25 35 2b 01 00 	bnd jmp *0x12b35(%rip)        # 1c940 <MPI_Recv_init@Base>
    9e0b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e10 <MPI_Cart_create@plt>:
    9e10:	f3 0f 1e fa          	endbr64 
    9e14:	f2 ff 25 2d 2b 01 00 	bnd jmp *0x12b2d(%rip)        # 1c948 <MPI_Cart_create@Base>
    9e1b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e20 <MPI_Comm_get_name@plt>:
    9e20:	f3 0f 1e fa          	endbr64 
    9e24:	f2 ff 25 25 2b 01 00 	bnd jmp *0x12b25(%rip)        # 1c950 <MPI_Comm_get_name@Base>
    9e2b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e30 <MPI_Info_dup@plt>:
    9e30:	f3 0f 1e fa          	endbr64 
    9e34:	f2 ff 25 1d 2b 01 00 	bnd jmp *0x12b1d(%rip)        # 1c958 <MPI_Info_dup@Base>
    9e3b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e40 <MPI_Irsend@plt>:
    9e40:	f3 0f 1e fa          	endbr64 
    9e44:	f2 ff 25 15 2b 01 00 	bnd jmp *0x12b15(%rip)        # 1c960 <MPI_Irsend@Base>
    9e4b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e50 <MPI_Bcast@plt>:
    9e50:	f3 0f 1e fa          	endbr64 
    9e54:	f2 ff 25 0d 2b 01 00 	bnd jmp *0x12b0d(%rip)        # 1c968 <MPI_Bcast@Base>
    9e5b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e60 <_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE@plt>:
    9e60:	f3 0f 1e fa          	endbr64 
    9e64:	f2 ff 25 05 2b 01 00 	bnd jmp *0x12b05(%rip)        # 1c970 <_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE@GLIBCXX_3.4.21>
    9e6b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e70 <MPI_Ssend@plt>:
    9e70:	f3 0f 1e fa          	endbr64 
    9e74:	f2 ff 25 fd 2a 01 00 	bnd jmp *0x12afd(%rip)        # 1c978 <MPI_Ssend@Base>
    9e7b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e80 <MPI_Group_range_incl@plt>:
    9e80:	f3 0f 1e fa          	endbr64 
    9e84:	f2 ff 25 f5 2a 01 00 	bnd jmp *0x12af5(%rip)        # 1c980 <MPI_Group_range_incl@Base>
    9e8b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009e90 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEE5closeEv@plt>:
    9e90:	f3 0f 1e fa          	endbr64 
    9e94:	f2 ff 25 ed 2a 01 00 	bnd jmp *0x12aed(%rip)        # 1c988 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEE5closeEv@GLIBCXX_3.4>
    9e9b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009ea0 <MPI_Op_create@plt>:
    9ea0:	f3 0f 1e fa          	endbr64 
    9ea4:	f2 ff 25 e5 2a 01 00 	bnd jmp *0x12ae5(%rip)        # 1c990 <MPI_Op_create@Base>
    9eab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009eb0 <MPI_Group_size@plt>:
    9eb0:	f3 0f 1e fa          	endbr64 
    9eb4:	f2 ff 25 dd 2a 01 00 	bnd jmp *0x12add(%rip)        # 1c998 <MPI_Group_size@Base>
    9ebb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009ec0 <MPI_Comm_disconnect@plt>:
    9ec0:	f3 0f 1e fa          	endbr64 
    9ec4:	f2 ff 25 d5 2a 01 00 	bnd jmp *0x12ad5(%rip)        # 1c9a0 <MPI_Comm_disconnect@Base>
    9ecb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009ed0 <_ZNSolsEm@plt>:
    9ed0:	f3 0f 1e fa          	endbr64 
    9ed4:	f2 ff 25 cd 2a 01 00 	bnd jmp *0x12acd(%rip)        # 1c9a8 <_ZNSolsEm@GLIBCXX_3.4>
    9edb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009ee0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmPKc@plt>:
    9ee0:	f3 0f 1e fa          	endbr64 
    9ee4:	f2 ff 25 c5 2a 01 00 	bnd jmp *0x12ac5(%rip)        # 1c9b0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmPKc@GLIBCXX_3.4.21>
    9eeb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009ef0 <MPI_Win_get_name@plt>:
    9ef0:	f3 0f 1e fa          	endbr64 
    9ef4:	f2 ff 25 bd 2a 01 00 	bnd jmp *0x12abd(%rip)        # 1c9b8 <MPI_Win_get_name@Base>
    9efb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f00 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEED1Ev@plt>:
    9f00:	f3 0f 1e fa          	endbr64 
    9f04:	f2 ff 25 b5 2a 01 00 	bnd jmp *0x12ab5(%rip)        # 1c9c0 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEED1Ev@GLIBCXX_3.4>
    9f0b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f10 <MPI_Pack_size@plt>:
    9f10:	f3 0f 1e fa          	endbr64 
    9f14:	f2 ff 25 ad 2a 01 00 	bnd jmp *0x12aad(%rip)        # 1c9c8 <MPI_Pack_size@Base>
    9f1b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f20 <MPI_Request_get_status@plt>:
    9f20:	f3 0f 1e fa          	endbr64 
    9f24:	f2 ff 25 a5 2a 01 00 	bnd jmp *0x12aa5(%rip)        # 1c9d0 <MPI_Request_get_status@Base>
    9f2b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f30 <MPI_Graph_neighbors@plt>:
    9f30:	f3 0f 1e fa          	endbr64 
    9f34:	f2 ff 25 9d 2a 01 00 	bnd jmp *0x12a9d(%rip)        # 1c9d8 <MPI_Graph_neighbors@Base>
    9f3b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f40 <MPI_Type_create_hvector@plt>:
    9f40:	f3 0f 1e fa          	endbr64 
    9f44:	f2 ff 25 95 2a 01 00 	bnd jmp *0x12a95(%rip)        # 1c9e0 <MPI_Type_create_hvector@Base>
    9f4b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f50 <MPI_Unpack_external@plt>:
    9f50:	f3 0f 1e fa          	endbr64 
    9f54:	f2 ff 25 8d 2a 01 00 	bnd jmp *0x12a8d(%rip)        # 1c9e8 <MPI_Unpack_external@Base>
    9f5b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f60 <MPI_Probe@plt>:
    9f60:	f3 0f 1e fa          	endbr64 
    9f64:	f2 ff 25 85 2a 01 00 	bnd jmp *0x12a85(%rip)        # 1c9f0 <MPI_Probe@Base>
    9f6b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f70 <_ZNSolsEPFRSt8ios_baseS0_E@plt>:
    9f70:	f3 0f 1e fa          	endbr64 
    9f74:	f2 ff 25 7d 2a 01 00 	bnd jmp *0x12a7d(%rip)        # 1c9f8 <_ZNSolsEPFRSt8ios_baseS0_E@GLIBCXX_3.4>
    9f7b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f80 <MPI_Send@plt>:
    9f80:	f3 0f 1e fa          	endbr64 
    9f84:	f2 ff 25 75 2a 01 00 	bnd jmp *0x12a75(%rip)        # 1ca00 <MPI_Send@Base>
    9f8b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009f90 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE13get_allocatorEv@plt>:
    9f90:	f3 0f 1e fa          	endbr64 
    9f94:	f2 ff 25 6d 2a 01 00 	bnd jmp *0x12a6d(%rip)        # 1ca08 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE13get_allocatorEv@GLIBCXX_3.4.21>
    9f9b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009fa0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE5c_strEv@plt>:
    9fa0:	f3 0f 1e fa          	endbr64 
    9fa4:	f2 ff 25 65 2a 01 00 	bnd jmp *0x12a65(%rip)        # 1ca10 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE5c_strEv@GLIBCXX_3.4.21>
    9fab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009fb0 <MPI_Info_set@plt>:
    9fb0:	f3 0f 1e fa          	endbr64 
    9fb4:	f2 ff 25 5d 2a 01 00 	bnd jmp *0x12a5d(%rip)        # 1ca18 <MPI_Info_set@Base>
    9fbb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009fc0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>:
    9fc0:	f3 0f 1e fa          	endbr64 
    9fc4:	f2 ff 25 55 2a 01 00 	bnd jmp *0x12a55(%rip)        # 1ca20 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@GLIBCXX_3.4.21>
    9fcb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009fd0 <MPI_Type_get_extent@plt>:
    9fd0:	f3 0f 1e fa          	endbr64 
    9fd4:	f2 ff 25 4d 2a 01 00 	bnd jmp *0x12a4d(%rip)        # 1ca28 <MPI_Type_get_extent@Base>
    9fdb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009fe0 <MPI_Win_post@plt>:
    9fe0:	f3 0f 1e fa          	endbr64 
    9fe4:	f2 ff 25 45 2a 01 00 	bnd jmp *0x12a45(%rip)        # 1ca30 <MPI_Win_post@Base>
    9feb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

0000000000009ff0 <MPI_Isend@plt>:
    9ff0:	f3 0f 1e fa          	endbr64 
    9ff4:	f2 ff 25 3d 2a 01 00 	bnd jmp *0x12a3d(%rip)        # 1ca38 <MPI_Isend@Base>
    9ffb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a000 <MPI_Unpack@plt>:
    a000:	f3 0f 1e fa          	endbr64 
    a004:	f2 ff 25 35 2a 01 00 	bnd jmp *0x12a35(%rip)        # 1ca40 <MPI_Unpack@Base>
    a00b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a010 <MPI_Iprobe@plt>:
    a010:	f3 0f 1e fa          	endbr64 
    a014:	f2 ff 25 2d 2a 01 00 	bnd jmp *0x12a2d(%rip)        # 1ca48 <MPI_Iprobe@Base>
    a01b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a020 <MPI_Comm_set_attr@plt>:
    a020:	f3 0f 1e fa          	endbr64 
    a024:	f2 ff 25 25 2a 01 00 	bnd jmp *0x12a25(%rip)        # 1ca50 <MPI_Comm_set_attr@Base>
    a02b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a030 <MPI_Type_get_envelope@plt>:
    a030:	f3 0f 1e fa          	endbr64 
    a034:	f2 ff 25 1d 2a 01 00 	bnd jmp *0x12a1d(%rip)        # 1ca58 <MPI_Type_get_envelope@Base>
    a03b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a040 <MPI_Alltoall@plt>:
    a040:	f3 0f 1e fa          	endbr64 
    a044:	f2 ff 25 15 2a 01 00 	bnd jmp *0x12a15(%rip)        # 1ca60 <MPI_Alltoall@Base>
    a04b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a050 <MPI_Cancel@plt>:
    a050:	f3 0f 1e fa          	endbr64 
    a054:	f2 ff 25 0d 2a 01 00 	bnd jmp *0x12a0d(%rip)        # 1ca68 <MPI_Cancel@Base>
    a05b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a060 <MPI_Initialized@plt>:
    a060:	f3 0f 1e fa          	endbr64 
    a064:	f2 ff 25 05 2a 01 00 	bnd jmp *0x12a05(%rip)        # 1ca70 <MPI_Initialized@Base>
    a06b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a070 <MPI_Scatterv@plt>:
    a070:	f3 0f 1e fa          	endbr64 
    a074:	f2 ff 25 fd 29 01 00 	bnd jmp *0x129fd(%rip)        # 1ca78 <MPI_Scatterv@Base>
    a07b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a080 <MPI_Issend@plt>:
    a080:	f3 0f 1e fa          	endbr64 
    a084:	f2 ff 25 f5 29 01 00 	bnd jmp *0x129f5(%rip)        # 1ca80 <MPI_Issend@Base>
    a08b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a090 <_ZN3MPI4CommC2Ev@plt>:
    a090:	f3 0f 1e fa          	endbr64 
    a094:	f2 ff 25 ed 29 01 00 	bnd jmp *0x129ed(%rip)        # 1ca88 <_ZN3MPI4CommC2Ev@Base>
    a09b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a0a0 <MPI_Graphdims_get@plt>:
    a0a0:	f3 0f 1e fa          	endbr64 
    a0a4:	f2 ff 25 e5 29 01 00 	bnd jmp *0x129e5(%rip)        # 1ca90 <MPI_Graphdims_get@Base>
    a0ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a0b0 <MPI_Info_get@plt>:
    a0b0:	f3 0f 1e fa          	endbr64 
    a0b4:	f2 ff 25 dd 29 01 00 	bnd jmp *0x129dd(%rip)        # 1ca98 <MPI_Info_get@Base>
    a0bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a0c0 <MPI_Cart_shift@plt>:
    a0c0:	f3 0f 1e fa          	endbr64 
    a0c4:	f2 ff 25 d5 29 01 00 	bnd jmp *0x129d5(%rip)        # 1caa0 <MPI_Cart_shift@Base>
    a0cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a0d0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@plt>:
    a0d0:	f3 0f 1e fa          	endbr64 
    a0d4:	f2 ff 25 cd 29 01 00 	bnd jmp *0x129cd(%rip)        # 1caa8 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@GLIBCXX_3.4.21>
    a0db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a0e0 <MPI_Type_get_true_extent@plt>:
    a0e0:	f3 0f 1e fa          	endbr64 
    a0e4:	f2 ff 25 c5 29 01 00 	bnd jmp *0x129c5(%rip)        # 1cab0 <MPI_Type_get_true_extent@Base>
    a0eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a0f0 <MPI_Get@plt>:
    a0f0:	f3 0f 1e fa          	endbr64 
    a0f4:	f2 ff 25 bd 29 01 00 	bnd jmp *0x129bd(%rip)        # 1cab8 <MPI_Get@Base>
    a0fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>:
    a100:	f3 0f 1e fa          	endbr64 
    a104:	f2 ff 25 b5 29 01 00 	bnd jmp *0x129b5(%rip)        # 1cac0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@GLIBCXX_3.4.21>
    a10b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a110 <MPI_Type_contiguous@plt>:
    a110:	f3 0f 1e fa          	endbr64 
    a114:	f2 ff 25 ad 29 01 00 	bnd jmp *0x129ad(%rip)        # 1cac8 <MPI_Type_contiguous@Base>
    a11b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a120 <MPI_Cart_sub@plt>:
    a120:	f3 0f 1e fa          	endbr64 
    a124:	f2 ff 25 a5 29 01 00 	bnd jmp *0x129a5(%rip)        # 1cad0 <MPI_Cart_sub@Base>
    a12b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a130 <__assert_fail@plt>:
    a130:	f3 0f 1e fa          	endbr64 
    a134:	f2 ff 25 9d 29 01 00 	bnd jmp *0x1299d(%rip)        # 1cad8 <__assert_fail@GLIBC_2.2.5>
    a13b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a140 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EmcRKS3_@plt>:
    a140:	f3 0f 1e fa          	endbr64 
    a144:	f2 ff 25 95 29 01 00 	bnd jmp *0x12995(%rip)        # 1cae0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EmcRKS3_@GLIBCXX_3.4.21>
    a14b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a150 <MPI_Scatter@plt>:
    a150:	f3 0f 1e fa          	endbr64 
    a154:	f2 ff 25 8d 29 01 00 	bnd jmp *0x1298d(%rip)        # 1cae8 <MPI_Scatter@Base>
    a15b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a160 <MPI_Win_get_group@plt>:
    a160:	f3 0f 1e fa          	endbr64 
    a164:	f2 ff 25 85 29 01 00 	bnd jmp *0x12985(%rip)        # 1caf0 <MPI_Win_get_group@Base>
    a16b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a170 <_ZNSo5flushEv@plt>:
    a170:	f3 0f 1e fa          	endbr64 
    a174:	f2 ff 25 7d 29 01 00 	bnd jmp *0x1297d(%rip)        # 1caf8 <_ZNSo5flushEv@GLIBCXX_3.4>
    a17b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a180 <MPI_Ssend_init@plt>:
    a180:	f3 0f 1e fa          	endbr64 
    a184:	f2 ff 25 75 29 01 00 	bnd jmp *0x12975(%rip)        # 1cb00 <MPI_Ssend_init@Base>
    a18b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a190 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE4sizeEv@plt>:
    a190:	f3 0f 1e fa          	endbr64 
    a194:	f2 ff 25 6d 29 01 00 	bnd jmp *0x1296d(%rip)        # 1cb08 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE4sizeEv@GLIBCXX_3.4.21>
    a19b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a1a0 <MPI_Comm_spawn@plt>:
    a1a0:	f3 0f 1e fa          	endbr64 
    a1a4:	f2 ff 25 65 29 01 00 	bnd jmp *0x12965(%rip)        # 1cb10 <MPI_Comm_spawn@Base>
    a1ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a1b0 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEC1Ev@plt>:
    a1b0:	f3 0f 1e fa          	endbr64 
    a1b4:	f2 ff 25 5d 29 01 00 	bnd jmp *0x1295d(%rip)        # 1cb18 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEC1Ev@GLIBCXX_3.4.26>
    a1bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a1c0 <MPI_Barrier@plt>:
    a1c0:	f3 0f 1e fa          	endbr64 
    a1c4:	f2 ff 25 55 29 01 00 	bnd jmp *0x12955(%rip)        # 1cb20 <MPI_Barrier@Base>
    a1cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a1d0 <MPI_Win_start@plt>:
    a1d0:	f3 0f 1e fa          	endbr64 
    a1d4:	f2 ff 25 4d 29 01 00 	bnd jmp *0x1294d(%rip)        # 1cb28 <MPI_Win_start@Base>
    a1db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a1e0 <MPI_Pack_external@plt>:
    a1e0:	f3 0f 1e fa          	endbr64 
    a1e4:	f2 ff 25 45 29 01 00 	bnd jmp *0x12945(%rip)        # 1cb30 <MPI_Pack_external@Base>
    a1eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a1f0 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEED1Ev@plt>:
    a1f0:	f3 0f 1e fa          	endbr64 
    a1f4:	f2 ff 25 3d 29 01 00 	bnd jmp *0x1293d(%rip)        # 1cb38 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEED1Ev@GLIBCXX_3.4.21>
    a1fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a200 <MPI_Rsend_init@plt>:
    a200:	f3 0f 1e fa          	endbr64 
    a204:	f2 ff 25 35 29 01 00 	bnd jmp *0x12935(%rip)        # 1cb40 <MPI_Rsend_init@Base>
    a20b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a210 <MPI_Group_incl@plt>:
    a210:	f3 0f 1e fa          	endbr64 
    a214:	f2 ff 25 2d 29 01 00 	bnd jmp *0x1292d(%rip)        # 1cb48 <MPI_Group_incl@Base>
    a21b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a220 <MPI_Bsend_init@plt>:
    a220:	f3 0f 1e fa          	endbr64 
    a224:	f2 ff 25 25 29 01 00 	bnd jmp *0x12925(%rip)        # 1cb50 <MPI_Bsend_init@Base>
    a22b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a230 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEE7is_openEv@plt>:
    a230:	f3 0f 1e fa          	endbr64 
    a234:	f2 ff 25 1d 29 01 00 	bnd jmp *0x1291d(%rip)        # 1cb58 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEE7is_openEv@GLIBCXX_3.4>
    a23b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a240 <MPI_Sendrecv@plt>:
    a240:	f3 0f 1e fa          	endbr64 
    a244:	f2 ff 25 15 29 01 00 	bnd jmp *0x12915(%rip)        # 1cb60 <MPI_Sendrecv@Base>
    a24b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a250 <MPI_Finalize@plt>:
    a250:	f3 0f 1e fa          	endbr64 
    a254:	f2 ff 25 0d 29 01 00 	bnd jmp *0x1290d(%rip)        # 1cb68 <MPI_Finalize@Base>
    a25b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a260 <sqrt@plt>:
    a260:	f3 0f 1e fa          	endbr64 
    a264:	f2 ff 25 05 29 01 00 	bnd jmp *0x12905(%rip)        # 1cb70 <sqrt@GLIBC_2.2.5>
    a26b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a270 <MPI_Request_free@plt>:
    a270:	f3 0f 1e fa          	endbr64 
    a274:	f2 ff 25 fd 28 01 00 	bnd jmp *0x128fd(%rip)        # 1cb78 <MPI_Request_free@Base>
    a27b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a280 <MPI_Gather@plt>:
    a280:	f3 0f 1e fa          	endbr64 
    a284:	f2 ff 25 f5 28 01 00 	bnd jmp *0x128f5(%rip)        # 1cb80 <MPI_Gather@Base>
    a28b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a290 <MPI_Grequest_complete@plt>:
    a290:	f3 0f 1e fa          	endbr64 
    a294:	f2 ff 25 ed 28 01 00 	bnd jmp *0x128ed(%rip)        # 1cb88 <MPI_Grequest_complete@Base>
    a29b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a2a0 <MPI_Type_create_darray@plt>:
    a2a0:	f3 0f 1e fa          	endbr64 
    a2a4:	f2 ff 25 e5 28 01 00 	bnd jmp *0x128e5(%rip)        # 1cb90 <MPI_Type_create_darray@Base>
    a2ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a2b0 <MPI_Win_wait@plt>:
    a2b0:	f3 0f 1e fa          	endbr64 
    a2b4:	f2 ff 25 dd 28 01 00 	bnd jmp *0x128dd(%rip)        # 1cb98 <MPI_Win_wait@Base>
    a2bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a2c0 <MPI_Rsend@plt>:
    a2c0:	f3 0f 1e fa          	endbr64 
    a2c4:	f2 ff 25 d5 28 01 00 	bnd jmp *0x128d5(%rip)        # 1cba0 <MPI_Rsend@Base>
    a2cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a2d0 <MPI_Recv@plt>:
    a2d0:	f3 0f 1e fa          	endbr64 
    a2d4:	f2 ff 25 cd 28 01 00 	bnd jmp *0x128cd(%rip)        # 1cba8 <MPI_Recv@Base>
    a2db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a2e0 <MPI_Graph_get@plt>:
    a2e0:	f3 0f 1e fa          	endbr64 
    a2e4:	f2 ff 25 c5 28 01 00 	bnd jmp *0x128c5(%rip)        # 1cbb0 <MPI_Graph_get@Base>
    a2eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a2f0 <MPI_Comm_group@plt>:
    a2f0:	f3 0f 1e fa          	endbr64 
    a2f4:	f2 ff 25 bd 28 01 00 	bnd jmp *0x128bd(%rip)        # 1cbb8 <MPI_Comm_group@Base>
    a2fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a300 <MPI_Comm_create@plt>:
    a300:	f3 0f 1e fa          	endbr64 
    a304:	f2 ff 25 b5 28 01 00 	bnd jmp *0x128b5(%rip)        # 1cbc0 <MPI_Comm_create@Base>
    a30b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a310 <MPI_Cart_get@plt>:
    a310:	f3 0f 1e fa          	endbr64 
    a314:	f2 ff 25 ad 28 01 00 	bnd jmp *0x128ad(%rip)        # 1cbc8 <MPI_Cart_get@Base>
    a31b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a320 <__cxa_atexit@plt>:
    a320:	f3 0f 1e fa          	endbr64 
    a324:	f2 ff 25 a5 28 01 00 	bnd jmp *0x128a5(%rip)        # 1cbd0 <__cxa_atexit@GLIBC_2.2.5>
    a32b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a330 <MPI_Type_get_attr@plt>:
    a330:	f3 0f 1e fa          	endbr64 
    a334:	f2 ff 25 9d 28 01 00 	bnd jmp *0x1289d(%rip)        # 1cbd8 <MPI_Type_get_attr@Base>
    a33b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a340 <MPI_Op_commutative@plt>:
    a340:	f3 0f 1e fa          	endbr64 
    a344:	f2 ff 25 95 28 01 00 	bnd jmp *0x12895(%rip)        # 1cbe0 <MPI_Op_commutative@Base>
    a34b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a350 <MPI_Graph_map@plt>:
    a350:	f3 0f 1e fa          	endbr64 
    a354:	f2 ff 25 8d 28 01 00 	bnd jmp *0x1288d(%rip)        # 1cbe8 <MPI_Graph_map@Base>
    a35b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a360 <MPI_Win_set_name@plt>:
    a360:	f3 0f 1e fa          	endbr64 
    a364:	f2 ff 25 85 28 01 00 	bnd jmp *0x12885(%rip)        # 1cbf0 <MPI_Win_set_name@Base>
    a36b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a370 <MPI_Test_cancelled@plt>:
    a370:	f3 0f 1e fa          	endbr64 
    a374:	f2 ff 25 7d 28 01 00 	bnd jmp *0x1287d(%rip)        # 1cbf8 <MPI_Test_cancelled@Base>
    a37b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a380 <MPI_Cart_map@plt>:
    a380:	f3 0f 1e fa          	endbr64 
    a384:	f2 ff 25 75 28 01 00 	bnd jmp *0x12875(%rip)        # 1cc00 <MPI_Cart_map@Base>
    a38b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a390 <MPI_Exscan@plt>:
    a390:	f3 0f 1e fa          	endbr64 
    a394:	f2 ff 25 6d 28 01 00 	bnd jmp *0x1286d(%rip)        # 1cc08 <MPI_Exscan@Base>
    a39b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a3a0 <MPI_Pack@plt>:
    a3a0:	f3 0f 1e fa          	endbr64 
    a3a4:	f2 ff 25 65 28 01 00 	bnd jmp *0x12865(%rip)        # 1cc10 <MPI_Pack@Base>
    a3ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a3b0 <__isoc99_sscanf@plt>:
    a3b0:	f3 0f 1e fa          	endbr64 
    a3b4:	f2 ff 25 5d 28 01 00 	bnd jmp *0x1285d(%rip)        # 1cc18 <__isoc99_sscanf@GLIBC_2.7>
    a3bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a3c0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>:
    a3c0:	f3 0f 1e fa          	endbr64 
    a3c4:	f2 ff 25 55 28 01 00 	bnd jmp *0x12855(%rip)        # 1cc20 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@GLIBCXX_3.4>
    a3cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a3d0 <_Znwm@plt>:
    a3d0:	f3 0f 1e fa          	endbr64 
    a3d4:	f2 ff 25 4d 28 01 00 	bnd jmp *0x1284d(%rip)        # 1cc28 <_Znwm@GLIBCXX_3.4>
    a3db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a3e0 <MPI_Errhandler_free@plt>:
    a3e0:	f3 0f 1e fa          	endbr64 
    a3e4:	f2 ff 25 45 28 01 00 	bnd jmp *0x12845(%rip)        # 1cc30 <MPI_Errhandler_free@Base>
    a3eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a3f0 <_ZdlPvm@plt>:
    a3f0:	f3 0f 1e fa          	endbr64 
    a3f4:	f2 ff 25 3d 28 01 00 	bnd jmp *0x1283d(%rip)        # 1cc38 <_ZdlPvm@CXXABI_1.3.9>
    a3fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a400 <_ZNSolsEPFRSoS_E@plt>:
    a400:	f3 0f 1e fa          	endbr64 
    a404:	f2 ff 25 35 28 01 00 	bnd jmp *0x12835(%rip)        # 1cc40 <_ZNSolsEPFRSoS_E@GLIBCXX_3.4>
    a40b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a410 <MPI_Abort@plt>:
    a410:	f3 0f 1e fa          	endbr64 
    a414:	f2 ff 25 2d 28 01 00 	bnd jmp *0x1282d(%rip)        # 1cc48 <MPI_Abort@Base>
    a41b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a420 <MPI_Intercomm_create@plt>:
    a420:	f3 0f 1e fa          	endbr64 
    a424:	f2 ff 25 25 28 01 00 	bnd jmp *0x12825(%rip)        # 1cc50 <MPI_Intercomm_create@Base>
    a42b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a430 <MPI_Cartdim_get@plt>:
    a430:	f3 0f 1e fa          	endbr64 
    a434:	f2 ff 25 1d 28 01 00 	bnd jmp *0x1281d(%rip)        # 1cc58 <MPI_Cartdim_get@Base>
    a43b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a440 <MPI_Cart_coords@plt>:
    a440:	f3 0f 1e fa          	endbr64 
    a444:	f2 ff 25 15 28 01 00 	bnd jmp *0x12815(%rip)        # 1cc60 <MPI_Cart_coords@Base>
    a44b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a450 <ompi_op_set_cxx_callback@plt>:
    a450:	f3 0f 1e fa          	endbr64 
    a454:	f2 ff 25 0d 28 01 00 	bnd jmp *0x1280d(%rip)        # 1cc68 <ompi_op_set_cxx_callback@Base>
    a45b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a460 <_ZNSaIcED1Ev@plt>:
    a460:	f3 0f 1e fa          	endbr64 
    a464:	f2 ff 25 05 28 01 00 	bnd jmp *0x12805(%rip)        # 1cc70 <_ZNSaIcED1Ev@GLIBCXX_3.4>
    a46b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a470 <MPI_Get_elements@plt>:
    a470:	f3 0f 1e fa          	endbr64 
    a474:	f2 ff 25 fd 27 01 00 	bnd jmp *0x127fd(%rip)        # 1cc78 <MPI_Get_elements@Base>
    a47b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a480 <MPI_Comm_delete_attr@plt>:
    a480:	f3 0f 1e fa          	endbr64 
    a484:	f2 ff 25 f5 27 01 00 	bnd jmp *0x127f5(%rip)        # 1cc80 <MPI_Comm_delete_attr@Base>
    a48b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a490 <MPI_Send_init@plt>:
    a490:	f3 0f 1e fa          	endbr64 
    a494:	f2 ff 25 ed 27 01 00 	bnd jmp *0x127ed(%rip)        # 1cc88 <MPI_Send_init@Base>
    a49b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a4a0 <MPI_Info_get_nkeys@plt>:
    a4a0:	f3 0f 1e fa          	endbr64 
    a4a4:	f2 ff 25 e5 27 01 00 	bnd jmp *0x127e5(%rip)        # 1cc90 <MPI_Info_get_nkeys@Base>
    a4ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a4b0 <MPI_Bsend@plt>:
    a4b0:	f3 0f 1e fa          	endbr64 
    a4b4:	f2 ff 25 dd 27 01 00 	bnd jmp *0x127dd(%rip)        # 1cc98 <MPI_Bsend@Base>
    a4bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a4c0 <MPI_Scan@plt>:
    a4c0:	f3 0f 1e fa          	endbr64 
    a4c4:	f2 ff 25 d5 27 01 00 	bnd jmp *0x127d5(%rip)        # 1cca0 <MPI_Scan@Base>
    a4cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a4d0 <__stack_chk_fail@plt>:
    a4d0:	f3 0f 1e fa          	endbr64 
    a4d4:	f2 ff 25 cd 27 01 00 	bnd jmp *0x127cd(%rip)        # 1cca8 <__stack_chk_fail@GLIBC_2.4>
    a4db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a4e0 <MPI_Comm_size@plt>:
    a4e0:	f3 0f 1e fa          	endbr64 
    a4e4:	f2 ff 25 c5 27 01 00 	bnd jmp *0x127c5(%rip)        # 1ccb0 <MPI_Comm_size@Base>
    a4eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a4f0 <free@plt>:
    a4f0:	f3 0f 1e fa          	endbr64 
    a4f4:	f2 ff 25 bd 27 01 00 	bnd jmp *0x127bd(%rip)        # 1ccb8 <free@GLIBC_2.2.5>
    a4fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a500 <MPI_Type_create_hindexed@plt>:
    a500:	f3 0f 1e fa          	endbr64 
    a504:	f2 ff 25 b5 27 01 00 	bnd jmp *0x127b5(%rip)        # 1ccc0 <MPI_Type_create_hindexed@Base>
    a50b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a510 <MPI_Reduce_local@plt>:
    a510:	f3 0f 1e fa          	endbr64 
    a514:	f2 ff 25 ad 27 01 00 	bnd jmp *0x127ad(%rip)        # 1ccc8 <MPI_Reduce_local@Base>
    a51b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a520 <MPI_Reduce_scatter@plt>:
    a520:	f3 0f 1e fa          	endbr64 
    a524:	f2 ff 25 a5 27 01 00 	bnd jmp *0x127a5(%rip)        # 1ccd0 <MPI_Reduce_scatter@Base>
    a52b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a530 <exit@plt>:
    a530:	f3 0f 1e fa          	endbr64 
    a534:	f2 ff 25 9d 27 01 00 	bnd jmp *0x1279d(%rip)        # 1ccd8 <exit@GLIBC_2.2.5>
    a53b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a540 <MPI_Win_lock@plt>:
    a540:	f3 0f 1e fa          	endbr64 
    a544:	f2 ff 25 95 27 01 00 	bnd jmp *0x12795(%rip)        # 1cce0 <MPI_Win_lock@Base>
    a54b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a550 <MPI_Allreduce@plt>:
    a550:	f3 0f 1e fa          	endbr64 
    a554:	f2 ff 25 8d 27 01 00 	bnd jmp *0x1278d(%rip)        # 1cce8 <MPI_Allreduce@Base>
    a55b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a560 <MPI_Allgatherv@plt>:
    a560:	f3 0f 1e fa          	endbr64 
    a564:	f2 ff 25 85 27 01 00 	bnd jmp *0x12785(%rip)        # 1ccf0 <MPI_Allgatherv@Base>
    a56b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a570 <atof@plt>:
    a570:	f3 0f 1e fa          	endbr64 
    a574:	f2 ff 25 7d 27 01 00 	bnd jmp *0x1277d(%rip)        # 1ccf8 <atof@GLIBC_2.2.5>
    a57b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a580 <MPI_Group_range_excl@plt>:
    a580:	f3 0f 1e fa          	endbr64 
    a584:	f2 ff 25 75 27 01 00 	bnd jmp *0x12775(%rip)        # 1cd00 <MPI_Group_range_excl@Base>
    a58b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a590 <_ZNKSt9basic_iosIcSt11char_traitsIcEEcvbEv@plt>:
    a590:	f3 0f 1e fa          	endbr64 
    a594:	f2 ff 25 6d 27 01 00 	bnd jmp *0x1276d(%rip)        # 1cd08 <_ZNKSt9basic_iosIcSt11char_traitsIcEEcvbEv@GLIBCXX_3.4.21>
    a59b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a5a0 <MPI_Type_commit@plt>:
    a5a0:	f3 0f 1e fa          	endbr64 
    a5a4:	f2 ff 25 65 27 01 00 	bnd jmp *0x12765(%rip)        # 1cd10 <MPI_Type_commit@Base>
    a5ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a5b0 <MPI_Put@plt>:
    a5b0:	f3 0f 1e fa          	endbr64 
    a5b4:	f2 ff 25 5d 27 01 00 	bnd jmp *0x1275d(%rip)        # 1cd18 <MPI_Put@Base>
    a5bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a5c0 <_ZdaPv@plt>:
    a5c0:	f3 0f 1e fa          	endbr64 
    a5c4:	f2 ff 25 55 27 01 00 	bnd jmp *0x12755(%rip)        # 1cd20 <_ZdaPv@GLIBCXX_3.4>
    a5cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a5d0 <MPI_Type_size@plt>:
    a5d0:	f3 0f 1e fa          	endbr64 
    a5d4:	f2 ff 25 4d 27 01 00 	bnd jmp *0x1274d(%rip)        # 1cd28 <MPI_Type_size@Base>
    a5db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a5e0 <MPI_Status_set_elements@plt>:
    a5e0:	f3 0f 1e fa          	endbr64 
    a5e4:	f2 ff 25 45 27 01 00 	bnd jmp *0x12745(%rip)        # 1cd30 <MPI_Status_set_elements@Base>
    a5eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a5f0 <MPI_Comm_get_errhandler@plt>:
    a5f0:	f3 0f 1e fa          	endbr64 
    a5f4:	f2 ff 25 3d 27 01 00 	bnd jmp *0x1273d(%rip)        # 1cd38 <MPI_Comm_get_errhandler@Base>
    a5fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a600 <MPI_Reduce@plt>:
    a600:	f3 0f 1e fa          	endbr64 
    a604:	f2 ff 25 35 27 01 00 	bnd jmp *0x12735(%rip)        # 1cd40 <MPI_Reduce@Base>
    a60b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a610 <MPI_Comm_remote_group@plt>:
    a610:	f3 0f 1e fa          	endbr64 
    a614:	f2 ff 25 2d 27 01 00 	bnd jmp *0x1272d(%rip)        # 1cd48 <MPI_Comm_remote_group@Base>
    a61b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a620 <malloc@plt>:
    a620:	f3 0f 1e fa          	endbr64 
    a624:	f2 ff 25 25 27 01 00 	bnd jmp *0x12725(%rip)        # 1cd50 <malloc@GLIBC_2.2.5>
    a62b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a630 <MPI_Wait@plt>:
    a630:	f3 0f 1e fa          	endbr64 
    a634:	f2 ff 25 1d 27 01 00 	bnd jmp *0x1271d(%rip)        # 1cd58 <MPI_Wait@Base>
    a63b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a640 <MPI_Comm_dup@plt>:
    a640:	f3 0f 1e fa          	endbr64 
    a644:	f2 ff 25 15 27 01 00 	bnd jmp *0x12715(%rip)        # 1cd60 <MPI_Comm_dup@Base>
    a64b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a650 <MPI_Sendrecv_replace@plt>:
    a650:	f3 0f 1e fa          	endbr64 
    a654:	f2 ff 25 0d 27 01 00 	bnd jmp *0x1270d(%rip)        # 1cd68 <MPI_Sendrecv_replace@Base>
    a65b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a660 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EPKcRKS3_@plt>:
    a660:	f3 0f 1e fa          	endbr64 
    a664:	f2 ff 25 05 27 01 00 	bnd jmp *0x12705(%rip)        # 1cd70 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EPKcRKS3_@GLIBCXX_3.4.21>
    a66b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a670 <MPI_Irecv@plt>:
    a670:	f3 0f 1e fa          	endbr64 
    a674:	f2 ff 25 fd 26 01 00 	bnd jmp *0x126fd(%rip)        # 1cd78 <MPI_Irecv@Base>
    a67b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a680 <MPI_Alltoallw@plt>:
    a680:	f3 0f 1e fa          	endbr64 
    a684:	f2 ff 25 f5 26 01 00 	bnd jmp *0x126f5(%rip)        # 1cd80 <MPI_Alltoallw@Base>
    a68b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a690 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1Ev@plt>:
    a690:	f3 0f 1e fa          	endbr64 
    a694:	f2 ff 25 ed 26 01 00 	bnd jmp *0x126ed(%rip)        # 1cd88 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1Ev@GLIBCXX_3.4.21>
    a69b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a6a0 <MPI_Type_get_name@plt>:
    a6a0:	f3 0f 1e fa          	endbr64 
    a6a4:	f2 ff 25 e5 26 01 00 	bnd jmp *0x126e5(%rip)        # 1cd90 <MPI_Type_get_name@Base>
    a6ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a6b0 <MPI_Type_create_indexed_block@plt>:
    a6b0:	f3 0f 1e fa          	endbr64 
    a6b4:	f2 ff 25 dd 26 01 00 	bnd jmp *0x126dd(%rip)        # 1cd98 <MPI_Type_create_indexed_block@Base>
    a6bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a6c0 <MPI_Type_get_contents@plt>:
    a6c0:	f3 0f 1e fa          	endbr64 
    a6c4:	f2 ff 25 d5 26 01 00 	bnd jmp *0x126d5(%rip)        # 1cda0 <MPI_Type_get_contents@Base>
    a6cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a6d0 <MPI_Comm_get_attr@plt>:
    a6d0:	f3 0f 1e fa          	endbr64 
    a6d4:	f2 ff 25 cd 26 01 00 	bnd jmp *0x126cd(%rip)        # 1cda8 <MPI_Comm_get_attr@Base>
    a6db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a6e0 <__cxa_throw_bad_array_new_length@plt>:
    a6e0:	f3 0f 1e fa          	endbr64 
    a6e4:	f2 ff 25 c5 26 01 00 	bnd jmp *0x126c5(%rip)        # 1cdb0 <__cxa_throw_bad_array_new_length@CXXABI_1.3.8>
    a6eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a6f0 <MPI_Comm_rank@plt>:
    a6f0:	f3 0f 1e fa          	endbr64 
    a6f4:	f2 ff 25 bd 26 01 00 	bnd jmp *0x126bd(%rip)        # 1cdb8 <MPI_Comm_rank@Base>
    a6fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a700 <MPI_Init@plt>:
    a700:	f3 0f 1e fa          	endbr64 
    a704:	f2 ff 25 b5 26 01 00 	bnd jmp *0x126b5(%rip)        # 1cdc0 <MPI_Init@Base>
    a70b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a710 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEE5closeEv@plt>:
    a710:	f3 0f 1e fa          	endbr64 
    a714:	f2 ff 25 ad 26 01 00 	bnd jmp *0x126ad(%rip)        # 1cdc8 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEE5closeEv@GLIBCXX_3.4>
    a71b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a720 <MPI_Info_delete@plt>:
    a720:	f3 0f 1e fa          	endbr64 
    a724:	f2 ff 25 a5 26 01 00 	bnd jmp *0x126a5(%rip)        # 1cdd0 <MPI_Info_delete@Base>
    a72b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a730 <MPI_Topo_test@plt>:
    a730:	f3 0f 1e fa          	endbr64 
    a734:	f2 ff 25 9d 26 01 00 	bnd jmp *0x1269d(%rip)        # 1cdd8 <MPI_Topo_test@Base>
    a73b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a740 <MPI_Get_count@plt>:
    a740:	f3 0f 1e fa          	endbr64 
    a744:	f2 ff 25 95 26 01 00 	bnd jmp *0x12695(%rip)        # 1cde0 <MPI_Get_count@Base>
    a74b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a750 <fprintf@plt>:
    a750:	f3 0f 1e fa          	endbr64 
    a754:	f2 ff 25 8d 26 01 00 	bnd jmp *0x1268d(%rip)        # 1cde8 <fprintf@GLIBC_2.2.5>
    a75b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a760 <MPI_Info_get_nthkey@plt>:
    a760:	f3 0f 1e fa          	endbr64 
    a764:	f2 ff 25 85 26 01 00 	bnd jmp *0x12685(%rip)        # 1cdf0 <MPI_Info_get_nthkey@Base>
    a76b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a770 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEED1Ev@plt>:
    a770:	f3 0f 1e fa          	endbr64 
    a774:	f2 ff 25 7d 26 01 00 	bnd jmp *0x1267d(%rip)        # 1cdf8 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEED1Ev@GLIBCXX_3.4>
    a77b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a780 <MPI_Win_set_errhandler@plt>:
    a780:	f3 0f 1e fa          	endbr64 
    a784:	f2 ff 25 75 26 01 00 	bnd jmp *0x12675(%rip)        # 1ce00 <MPI_Win_set_errhandler@Base>
    a78b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a790 <MPI_Comm_remote_size@plt>:
    a790:	f3 0f 1e fa          	endbr64 
    a794:	f2 ff 25 6d 26 01 00 	bnd jmp *0x1266d(%rip)        # 1ce08 <MPI_Comm_remote_size@Base>
    a79b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a7a0 <MPI_Type_dup@plt>:
    a7a0:	f3 0f 1e fa          	endbr64 
    a7a4:	f2 ff 25 65 26 01 00 	bnd jmp *0x12665(%rip)        # 1ce10 <MPI_Type_dup@Base>
    a7ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a7b0 <_ZNSt8ios_base4InitC1Ev@plt>:
    a7b0:	f3 0f 1e fa          	endbr64 
    a7b4:	f2 ff 25 5d 26 01 00 	bnd jmp *0x1265d(%rip)        # 1ce18 <_ZNSt8ios_base4InitC1Ev@GLIBCXX_3.4>
    a7bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a7c0 <puts@plt>:
    a7c0:	f3 0f 1e fa          	endbr64 
    a7c4:	f2 ff 25 55 26 01 00 	bnd jmp *0x12655(%rip)        # 1ce20 <puts@GLIBC_2.2.5>
    a7cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a7d0 <MPI_Comm_set_errhandler@plt>:
    a7d0:	f3 0f 1e fa          	endbr64 
    a7d4:	f2 ff 25 4d 26 01 00 	bnd jmp *0x1264d(%rip)        # 1ce28 <MPI_Comm_set_errhandler@Base>
    a7db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a7e0 <MPI_Ibsend@plt>:
    a7e0:	f3 0f 1e fa          	endbr64 
    a7e4:	f2 ff 25 45 26 01 00 	bnd jmp *0x12645(%rip)        # 1ce30 <MPI_Ibsend@Base>
    a7eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a7f0 <MPI_Alltoallv@plt>:
    a7f0:	f3 0f 1e fa          	endbr64 
    a7f4:	f2 ff 25 3d 26 01 00 	bnd jmp *0x1263d(%rip)        # 1ce38 <MPI_Alltoallv@Base>
    a7fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a800 <MPI_Op_free@plt>:
    a800:	f3 0f 1e fa          	endbr64 
    a804:	f2 ff 25 35 26 01 00 	bnd jmp *0x12635(%rip)        # 1ce40 <MPI_Op_free@Base>
    a80b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a810 <strtok@plt>:
    a810:	f3 0f 1e fa          	endbr64 
    a814:	f2 ff 25 2d 26 01 00 	bnd jmp *0x1262d(%rip)        # 1ce48 <strtok@GLIBC_2.2.5>
    a81b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a820 <MPI_Comm_connect@plt>:
    a820:	f3 0f 1e fa          	endbr64 
    a824:	f2 ff 25 25 26 01 00 	bnd jmp *0x12625(%rip)        # 1ce50 <MPI_Comm_connect@Base>
    a82b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a830 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmRKS4_@plt>:
    a830:	f3 0f 1e fa          	endbr64 
    a834:	f2 ff 25 1d 26 01 00 	bnd jmp *0x1261d(%rip)        # 1ce58 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmRKS4_@GLIBCXX_3.4.21>
    a83b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a840 <MPI_Comm_free@plt>:
    a840:	f3 0f 1e fa          	endbr64 
    a844:	f2 ff 25 15 26 01 00 	bnd jmp *0x12615(%rip)        # 1ce60 <MPI_Comm_free@Base>
    a84b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a850 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEE7is_openEv@plt>:
    a850:	f3 0f 1e fa          	endbr64 
    a854:	f2 ff 25 0d 26 01 00 	bnd jmp *0x1260d(%rip)        # 1ce68 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEE7is_openEv@GLIBCXX_3.4>
    a85b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a860 <MPI_Comm_split@plt>:
    a860:	f3 0f 1e fa          	endbr64 
    a864:	f2 ff 25 05 26 01 00 	bnd jmp *0x12605(%rip)        # 1ce70 <MPI_Comm_split@Base>
    a86b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a870 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEEC1ERKNSt7__cxx1112basic_stringIcS1_SaIcEEESt13_Ios_Openmode@plt>:
    a870:	f3 0f 1e fa          	endbr64 
    a874:	f2 ff 25 fd 25 01 00 	bnd jmp *0x125fd(%rip)        # 1ce78 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEEC1ERKNSt7__cxx1112basic_stringIcS1_SaIcEEESt13_Ios_Openmode@GLIBCXX_3.4.21>
    a87b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a880 <MPI_Type_create_subarray@plt>:
    a880:	f3 0f 1e fa          	endbr64 
    a884:	f2 ff 25 f5 25 01 00 	bnd jmp *0x125f5(%rip)        # 1ce80 <MPI_Type_create_subarray@Base>
    a88b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a890 <MPI_Comm_spawn_multiple@plt>:
    a890:	f3 0f 1e fa          	endbr64 
    a894:	f2 ff 25 ed 25 01 00 	bnd jmp *0x125ed(%rip)        # 1ce88 <MPI_Comm_spawn_multiple@Base>
    a89b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a8a0 <realloc@plt>:
    a8a0:	f3 0f 1e fa          	endbr64 
    a8a4:	f2 ff 25 e5 25 01 00 	bnd jmp *0x125e5(%rip)        # 1ce90 <realloc@GLIBC_2.2.5>
    a8ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a8b0 <MPI_Test@plt>:
    a8b0:	f3 0f 1e fa          	endbr64 
    a8b4:	f2 ff 25 dd 25 01 00 	bnd jmp *0x125dd(%rip)        # 1ce98 <MPI_Test@Base>
    a8bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a8c0 <MPI_Info_get_valuelen@plt>:
    a8c0:	f3 0f 1e fa          	endbr64 
    a8c4:	f2 ff 25 d5 25 01 00 	bnd jmp *0x125d5(%rip)        # 1cea0 <MPI_Info_get_valuelen@Base>
    a8cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a8d0 <MPI_Win_test@plt>:
    a8d0:	f3 0f 1e fa          	endbr64 
    a8d4:	f2 ff 25 cd 25 01 00 	bnd jmp *0x125cd(%rip)        # 1cea8 <MPI_Win_test@Base>
    a8db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a8e0 <MPI_Type_vector@plt>:
    a8e0:	f3 0f 1e fa          	endbr64 
    a8e4:	f2 ff 25 c5 25 01 00 	bnd jmp *0x125c5(%rip)        # 1ceb0 <MPI_Type_vector@Base>
    a8eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a8f0 <MPI_Start@plt>:
    a8f0:	f3 0f 1e fa          	endbr64 
    a8f4:	f2 ff 25 bd 25 01 00 	bnd jmp *0x125bd(%rip)        # 1ceb8 <MPI_Start@Base>
    a8fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a900 <MPI_Win_call_errhandler@plt>:
    a900:	f3 0f 1e fa          	endbr64 
    a904:	f2 ff 25 b5 25 01 00 	bnd jmp *0x125b5(%rip)        # 1cec0 <MPI_Win_call_errhandler@Base>
    a90b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a910 <MPI_Type_indexed@plt>:
    a910:	f3 0f 1e fa          	endbr64 
    a914:	f2 ff 25 ad 25 01 00 	bnd jmp *0x125ad(%rip)        # 1cec8 <MPI_Type_indexed@Base>
    a91b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a920 <MPI_Comm_set_name@plt>:
    a920:	f3 0f 1e fa          	endbr64 
    a924:	f2 ff 25 a5 25 01 00 	bnd jmp *0x125a5(%rip)        # 1ced0 <MPI_Comm_set_name@Base>
    a92b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a930 <MPI_Win_fence@plt>:
    a930:	f3 0f 1e fa          	endbr64 
    a934:	f2 ff 25 9d 25 01 00 	bnd jmp *0x1259d(%rip)        # 1ced8 <MPI_Win_fence@Base>
    a93b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a940 <MPI_Win_complete@plt>:
    a940:	f3 0f 1e fa          	endbr64 
    a944:	f2 ff 25 95 25 01 00 	bnd jmp *0x12595(%rip)        # 1cee0 <MPI_Win_complete@Base>
    a94b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a950 <_ZNSolsEi@plt>:
    a950:	f3 0f 1e fa          	endbr64 
    a954:	f2 ff 25 8d 25 01 00 	bnd jmp *0x1258d(%rip)        # 1cee8 <_ZNSolsEi@GLIBCXX_3.4>
    a95b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a960 <MPI_Group_free@plt>:
    a960:	f3 0f 1e fa          	endbr64 
    a964:	f2 ff 25 85 25 01 00 	bnd jmp *0x12585(%rip)        # 1cef0 <MPI_Group_free@Base>
    a96b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a970 <MPI_Type_create_resized@plt>:
    a970:	f3 0f 1e fa          	endbr64 
    a974:	f2 ff 25 7d 25 01 00 	bnd jmp *0x1257d(%rip)        # 1cef8 <MPI_Type_create_resized@Base>
    a97b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a980 <_Unwind_Resume@plt>:
    a980:	f3 0f 1e fa          	endbr64 
    a984:	f2 ff 25 75 25 01 00 	bnd jmp *0x12575(%rip)        # 1cf00 <_Unwind_Resume@GCC_3.0>
    a98b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a990 <_ZNSaIcEC1Ev@plt>:
    a990:	f3 0f 1e fa          	endbr64 
    a994:	f2 ff 25 6d 25 01 00 	bnd jmp *0x1256d(%rip)        # 1cf08 <_ZNSaIcEC1Ev@GLIBCXX_3.4>
    a99b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a9a0 <MPI_Cart_rank@plt>:
    a9a0:	f3 0f 1e fa          	endbr64 
    a9a4:	f2 ff 25 65 25 01 00 	bnd jmp *0x12565(%rip)        # 1cf10 <MPI_Cart_rank@Base>
    a9ab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a9b0 <MPI_Comm_accept@plt>:
    a9b0:	f3 0f 1e fa          	endbr64 
    a9b4:	f2 ff 25 5d 25 01 00 	bnd jmp *0x1255d(%rip)        # 1cf18 <MPI_Comm_accept@Base>
    a9bb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a9c0 <MPI_Win_set_attr@plt>:
    a9c0:	f3 0f 1e fa          	endbr64 
    a9c4:	f2 ff 25 55 25 01 00 	bnd jmp *0x12555(%rip)        # 1cf20 <MPI_Win_set_attr@Base>
    a9cb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a9d0 <MPI_Win_delete_attr@plt>:
    a9d0:	f3 0f 1e fa          	endbr64 
    a9d4:	f2 ff 25 4d 25 01 00 	bnd jmp *0x1254d(%rip)        # 1cf28 <MPI_Win_delete_attr@Base>
    a9db:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a9e0 <MPI_Type_delete_attr@plt>:
    a9e0:	f3 0f 1e fa          	endbr64 
    a9e4:	f2 ff 25 45 25 01 00 	bnd jmp *0x12545(%rip)        # 1cf30 <MPI_Type_delete_attr@Base>
    a9eb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000a9f0 <MPI_Gatherv@plt>:
    a9f0:	f3 0f 1e fa          	endbr64 
    a9f4:	f2 ff 25 3d 25 01 00 	bnd jmp *0x1253d(%rip)        # 1cf38 <MPI_Gatherv@Base>
    a9fb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa00 <MPI_Group_rank@plt>:
    aa00:	f3 0f 1e fa          	endbr64 
    aa04:	f2 ff 25 35 25 01 00 	bnd jmp *0x12535(%rip)        # 1cf40 <MPI_Group_rank@Base>
    aa0b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa10 <MPI_Intercomm_merge@plt>:
    aa10:	f3 0f 1e fa          	endbr64 
    aa14:	f2 ff 25 2d 25 01 00 	bnd jmp *0x1252d(%rip)        # 1cf48 <MPI_Intercomm_merge@Base>
    aa1b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa20 <MPI_Win_get_errhandler@plt>:
    aa20:	f3 0f 1e fa          	endbr64 
    aa24:	f2 ff 25 25 25 01 00 	bnd jmp *0x12525(%rip)        # 1cf50 <MPI_Win_get_errhandler@Base>
    aa2b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa30 <MPI_Accumulate@plt>:
    aa30:	f3 0f 1e fa          	endbr64 
    aa34:	f2 ff 25 1d 25 01 00 	bnd jmp *0x1251d(%rip)        # 1cf58 <MPI_Accumulate@Base>
    aa3b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa40 <_ZNSolsEd@plt>:
    aa40:	f3 0f 1e fa          	endbr64 
    aa44:	f2 ff 25 15 25 01 00 	bnd jmp *0x12515(%rip)        # 1cf60 <_ZNSolsEd@GLIBCXX_3.4>
    aa4b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa50 <MPI_Type_set_attr@plt>:
    aa50:	f3 0f 1e fa          	endbr64 
    aa54:	f2 ff 25 0d 25 01 00 	bnd jmp *0x1250d(%rip)        # 1cf68 <MPI_Type_set_attr@Base>
    aa5b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa60 <MPI_Status_set_cancelled@plt>:
    aa60:	f3 0f 1e fa          	endbr64 
    aa64:	f2 ff 25 05 25 01 00 	bnd jmp *0x12505(%rip)        # 1cf70 <MPI_Status_set_cancelled@Base>
    aa6b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa70 <MPI_Graph_neighbors_count@plt>:
    aa70:	f3 0f 1e fa          	endbr64 
    aa74:	f2 ff 25 fd 24 01 00 	bnd jmp *0x124fd(%rip)        # 1cf78 <MPI_Graph_neighbors_count@Base>
    aa7b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa80 <MPI_Info_free@plt>:
    aa80:	f3 0f 1e fa          	endbr64 
    aa84:	f2 ff 25 f5 24 01 00 	bnd jmp *0x124f5(%rip)        # 1cf80 <MPI_Info_free@Base>
    aa8b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aa90 <MPI_Pack_external_size@plt>:
    aa90:	f3 0f 1e fa          	endbr64 
    aa94:	f2 ff 25 ed 24 01 00 	bnd jmp *0x124ed(%rip)        # 1cf88 <MPI_Pack_external_size@Base>
    aa9b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aaa0 <MPI_Type_set_name@plt>:
    aaa0:	f3 0f 1e fa          	endbr64 
    aaa4:	f2 ff 25 e5 24 01 00 	bnd jmp *0x124e5(%rip)        # 1cf90 <MPI_Type_set_name@Base>
    aaab:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aab0 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEEC1ERKNSt7__cxx1112basic_stringIcS1_SaIcEEESt13_Ios_Openmode@plt>:
    aab0:	f3 0f 1e fa          	endbr64 
    aab4:	f2 ff 25 dd 24 01 00 	bnd jmp *0x124dd(%rip)        # 1cf98 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEEC1ERKNSt7__cxx1112basic_stringIcS1_SaIcEEESt13_Ios_Openmode@GLIBCXX_3.4.21>
    aabb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aac0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendERKS4_@plt>:
    aac0:	f3 0f 1e fa          	endbr64 
    aac4:	f2 ff 25 d5 24 01 00 	bnd jmp *0x124d5(%rip)        # 1cfa0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendERKS4_@GLIBCXX_3.4.21>
    aacb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aad0 <MPI_Graph_create@plt>:
    aad0:	f3 0f 1e fa          	endbr64 
    aad4:	f2 ff 25 cd 24 01 00 	bnd jmp *0x124cd(%rip)        # 1cfa8 <MPI_Graph_create@Base>
    aadb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aae0 <MPI_Comm_test_inter@plt>:
    aae0:	f3 0f 1e fa          	endbr64 
    aae4:	f2 ff 25 c5 24 01 00 	bnd jmp *0x124c5(%rip)        # 1cfb0 <MPI_Comm_test_inter@Base>
    aaeb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

000000000000aaf0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEixEm@plt>:
    aaf0:	f3 0f 1e fa          	endbr64 
    aaf4:	f2 ff 25 bd 24 01 00 	bnd jmp *0x124bd(%rip)        # 1cfb8 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEixEm@GLIBCXX_3.4.21>
    aafb:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)

Disassembly of section .text:

000000000000ab00 <_start>:
    ab00:	f3 0f 1e fa          	endbr64 
    ab04:	31 ed                	xor    %ebp,%ebp
    ab06:	49 89 d1             	mov    %rdx,%r9
    ab09:	5e                   	pop    %rsi
    ab0a:	48 89 e2             	mov    %rsp,%rdx
    ab0d:	48 83 e4 f0          	and    $0xfffffffffffffff0,%rsp
    ab11:	50                   	push   %rax
    ab12:	54                   	push   %rsp
    ab13:	45 31 c0             	xor    %r8d,%r8d
    ab16:	31 c9                	xor    %ecx,%ecx
    ab18:	48 8d 3d ca 00 00 00 	lea    0xca(%rip),%rdi        # abe9 <main>
    ab1f:	ff 15 b3 24 01 00    	call   *0x124b3(%rip)        # 1cfd8 <__libc_start_main@GLIBC_2.34>
    ab25:	f4                   	hlt    
    ab26:	66 2e 0f 1f 84 00 00 	cs nopw 0x0(%rax,%rax,1)
    ab2d:	00 00 00 

000000000000ab30 <deregister_tm_clones>:
    ab30:	48 8d 3d e1 24 01 00 	lea    0x124e1(%rip),%rdi        # 1d018 <__TMC_END__>
    ab37:	48 8d 05 da 24 01 00 	lea    0x124da(%rip),%rax        # 1d018 <__TMC_END__>
    ab3e:	48 39 f8             	cmp    %rdi,%rax
    ab41:	74 15                	je     ab58 <deregister_tm_clones+0x28>
    ab43:	48 8b 05 96 24 01 00 	mov    0x12496(%rip),%rax        # 1cfe0 <_ITM_deregisterTMCloneTable@Base>
    ab4a:	48 85 c0             	test   %rax,%rax
    ab4d:	74 09                	je     ab58 <deregister_tm_clones+0x28>
    ab4f:	ff e0                	jmp    *%rax
    ab51:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)
    ab58:	c3                   	ret    
    ab59:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)

000000000000ab60 <register_tm_clones>:
    ab60:	48 8d 3d b1 24 01 00 	lea    0x124b1(%rip),%rdi        # 1d018 <__TMC_END__>
    ab67:	48 8d 35 aa 24 01 00 	lea    0x124aa(%rip),%rsi        # 1d018 <__TMC_END__>
    ab6e:	48 29 fe             	sub    %rdi,%rsi
    ab71:	48 89 f0             	mov    %rsi,%rax
    ab74:	48 c1 ee 3f          	shr    $0x3f,%rsi
    ab78:	48 c1 f8 03          	sar    $0x3,%rax
    ab7c:	48 01 c6             	add    %rax,%rsi
    ab7f:	48 d1 fe             	sar    %rsi
    ab82:	74 14                	je     ab98 <register_tm_clones+0x38>
    ab84:	48 8b 05 65 24 01 00 	mov    0x12465(%rip),%rax        # 1cff0 <_ITM_registerTMCloneTable@Base>
    ab8b:	48 85 c0             	test   %rax,%rax
    ab8e:	74 08                	je     ab98 <register_tm_clones+0x38>
    ab90:	ff e0                	jmp    *%rax
    ab92:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)
    ab98:	c3                   	ret    
    ab99:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)

000000000000aba0 <__do_global_dtors_aux>:
    aba0:	f3 0f 1e fa          	endbr64 
    aba4:	80 3d 7d 36 01 00 00 	cmpb   $0x0,0x1367d(%rip)        # 1e228 <completed.0>
    abab:	75 2b                	jne    abd8 <__do_global_dtors_aux+0x38>
    abad:	55                   	push   %rbp
    abae:	48 83 3d 12 24 01 00 	cmpq   $0x0,0x12412(%rip)        # 1cfc8 <__cxa_finalize@GLIBC_2.2.5>
    abb5:	00 
    abb6:	48 89 e5             	mov    %rsp,%rbp
    abb9:	74 0c                	je     abc7 <__do_global_dtors_aux+0x27>
    abbb:	48 8b 3d 46 24 01 00 	mov    0x12446(%rip),%rdi        # 1d008 <__dso_handle>
    abc2:	e8 c9 f1 ff ff       	call   9d90 <__cxa_finalize@plt>
    abc7:	e8 64 ff ff ff       	call   ab30 <deregister_tm_clones>
    abcc:	c6 05 55 36 01 00 01 	movb   $0x1,0x13655(%rip)        # 1e228 <completed.0>
    abd3:	5d                   	pop    %rbp
    abd4:	c3                   	ret    
    abd5:	0f 1f 00             	nopl   (%rax)
    abd8:	c3                   	ret    
    abd9:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)

000000000000abe0 <frame_dummy>:
    abe0:	f3 0f 1e fa          	endbr64 
    abe4:	e9 77 ff ff ff       	jmp    ab60 <register_tm_clones>

000000000000abe9 <main>:
    abe9:	f3 0f 1e fa          	endbr64 
    abed:	55                   	push   %rbp
    abee:	48 89 e5             	mov    %rsp,%rbp
    abf1:	41 55                	push   %r13
    abf3:	41 54                	push   %r12
    abf5:	53                   	push   %rbx
    abf6:	48 81 ec 88 01 00 00 	sub    $0x188,%rsp
    abfd:	89 bd 6c fe ff ff    	mov    %edi,-0x194(%rbp)
    ac03:	48 89 b5 60 fe ff ff 	mov    %rsi,-0x1a0(%rbp)
    ac0a:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ac11:	00 00 
    ac13:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
    ac17:	31 c0                	xor    %eax,%eax
    ac19:	c7 85 90 fe ff ff 03 	movl   $0x3,-0x170(%rbp)
    ac20:	00 00 00 
    ac23:	c7 85 94 fe ff ff 03 	movl   $0x3,-0x16c(%rbp)
    ac2a:	00 00 00 
    ac2d:	48 8d 95 60 fe ff ff 	lea    -0x1a0(%rbp),%rdx
    ac34:	48 8d 85 6c fe ff ff 	lea    -0x194(%rbp),%rax
    ac3b:	48 89 d6             	mov    %rdx,%rsi
    ac3e:	48 89 c7             	mov    %rax,%rdi
    ac41:	e8 ba fa ff ff       	call   a700 <MPI_Init@plt>
    ac46:	48 8d 85 c0 fe ff ff 	lea    -0x140(%rbp),%rax
    ac4d:	ba 03 00 00 00       	mov    $0x3,%edx
    ac52:	be 03 00 00 00       	mov    $0x3,%esi
    ac57:	48 89 c7             	mov    %rax,%rdi
    ac5a:	e8 54 78 00 00       	call   124b3 <_Z18setupCommunicatorsP10topology_tii>
    ac5f:	48 8d 85 c0 fe ff ff 	lea    -0x140(%rbp),%rax
    ac66:	48 89 85 b0 fe ff ff 	mov    %rax,-0x150(%rbp)
    ac6d:	48 8d 85 c0 fe ff ff 	lea    -0x140(%rbp),%rax
    ac74:	48 83 c0 04          	add    $0x4,%rax
    ac78:	48 89 85 b8 fe ff ff 	mov    %rax,-0x148(%rbp)
    ac7f:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    ac86:	48 89 c7             	mov    %rax,%rdi
    ac89:	e8 02 fd ff ff       	call   a990 <_ZNSaIcEC1Ev@plt>
    ac8e:	48 8d 95 a8 fe ff ff 	lea    -0x158(%rbp),%rdx
    ac95:	48 8d 85 f0 fe ff ff 	lea    -0x110(%rbp),%rax
    ac9c:	48 8d 0d a5 b5 00 00 	lea    0xb5a5(%rip),%rcx        # 16248 <_ZN3MPIL11LOCK_SHAREDE+0x8>
    aca3:	48 89 ce             	mov    %rcx,%rsi
    aca6:	48 89 c7             	mov    %rax,%rdi
    aca9:	e8 b2 f9 ff ff       	call   a660 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EPKcRKS3_@plt>
    acae:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    acb5:	48 89 c7             	mov    %rax,%rdi
    acb8:	e8 a3 f7 ff ff       	call   a460 <_ZNSaIcED1Ev@plt>
    acbd:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    acc4:	48 89 c7             	mov    %rax,%rdi
    acc7:	e8 c4 fc ff ff       	call   a990 <_ZNSaIcEC1Ev@plt>
    accc:	48 8d 95 a8 fe ff ff 	lea    -0x158(%rbp),%rdx
    acd3:	48 8d 85 10 ff ff ff 	lea    -0xf0(%rbp),%rax
    acda:	48 8d 0d bf b5 00 00 	lea    0xb5bf(%rip),%rcx        # 162a0 <_ZN3MPIL11LOCK_SHAREDE+0x60>
    ace1:	48 89 ce             	mov    %rcx,%rsi
    ace4:	48 89 c7             	mov    %rax,%rdi
    ace7:	e8 74 f9 ff ff       	call   a660 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EPKcRKS3_@plt>
    acec:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    acf3:	48 89 c7             	mov    %rax,%rdi
    acf6:	e8 65 f7 ff ff       	call   a460 <_ZNSaIcED1Ev@plt>
    acfb:	48 8b 85 b8 fe ff ff 	mov    -0x148(%rbp),%rax
    ad02:	8b 10                	mov    (%rax),%edx
    ad04:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    ad08:	89 d6                	mov    %edx,%esi
    ad0a:	48 89 c7             	mov    %rax,%rdi
    ad0d:	e8 3c 09 00 00       	call   b64e <_ZNSt7__cxx119to_stringEi>
    ad12:	48 8b 85 b0 fe ff ff 	mov    -0x150(%rbp),%rax
    ad19:	8b 10                	mov    (%rax),%edx
    ad1b:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    ad22:	89 d6                	mov    %edx,%esi
    ad24:	48 89 c7             	mov    %rax,%rdi
    ad27:	e8 22 09 00 00       	call   b64e <_ZNSt7__cxx119to_stringEi>
    ad2c:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    ad30:	48 8d 95 70 ff ff ff 	lea    -0x90(%rbp),%rdx
    ad37:	48 8d 0d ba b5 00 00 	lea    0xb5ba(%rip),%rcx        # 162f8 <_ZN3MPIL11LOCK_SHAREDE+0xb8>
    ad3e:	48 89 ce             	mov    %rcx,%rsi
    ad41:	48 89 c7             	mov    %rax,%rdi
    ad44:	e8 3a 5b 00 00       	call   10883 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEPKS5_OS8_>
    ad49:	48 8d 85 30 ff ff ff 	lea    -0xd0(%rbp),%rax
    ad50:	48 8d 55 b0          	lea    -0x50(%rbp),%rdx
    ad54:	48 8d 4d 90          	lea    -0x70(%rbp),%rcx
    ad58:	48 89 ce             	mov    %rcx,%rsi
    ad5b:	48 89 c7             	mov    %rax,%rdi
    ad5e:	e8 a7 59 00 00       	call   1070a <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_>
    ad63:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    ad67:	48 89 c7             	mov    %rax,%rdi
    ad6a:	e8 91 f3 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    ad6f:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    ad76:	48 89 c7             	mov    %rax,%rdi
    ad79:	e8 82 f3 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    ad7e:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    ad82:	48 89 c7             	mov    %rax,%rdi
    ad85:	e8 76 f3 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    ad8a:	48 8b 85 b8 fe ff ff 	mov    -0x148(%rbp),%rax
    ad91:	8b 10                	mov    (%rax),%edx
    ad93:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    ad97:	89 d6                	mov    %edx,%esi
    ad99:	48 89 c7             	mov    %rax,%rdi
    ad9c:	e8 ad 08 00 00       	call   b64e <_ZNSt7__cxx119to_stringEi>
    ada1:	48 8b 85 b0 fe ff ff 	mov    -0x150(%rbp),%rax
    ada8:	8b 10                	mov    (%rax),%edx
    adaa:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    adb1:	89 d6                	mov    %edx,%esi
    adb3:	48 89 c7             	mov    %rax,%rdi
    adb6:	e8 93 08 00 00       	call   b64e <_ZNSt7__cxx119to_stringEi>
    adbb:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    adbf:	48 8d 95 70 ff ff ff 	lea    -0x90(%rbp),%rdx
    adc6:	48 8d 0d 83 b5 00 00 	lea    0xb583(%rip),%rcx        # 16350 <_ZN3MPIL11LOCK_SHAREDE+0x110>
    adcd:	48 89 ce             	mov    %rcx,%rsi
    add0:	48 89 c7             	mov    %rax,%rdi
    add3:	e8 ab 5a 00 00       	call   10883 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEPKS5_OS8_>
    add8:	48 8d 85 50 ff ff ff 	lea    -0xb0(%rbp),%rax
    addf:	48 8d 55 b0          	lea    -0x50(%rbp),%rdx
    ade3:	48 8d 4d 90          	lea    -0x70(%rbp),%rcx
    ade7:	48 89 ce             	mov    %rcx,%rsi
    adea:	48 89 c7             	mov    %rax,%rdi
    aded:	e8 18 59 00 00       	call   1070a <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_>
    adf2:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    adf6:	48 89 c7             	mov    %rax,%rdi
    adf9:	e8 02 f3 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    adfe:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    ae05:	48 89 c7             	mov    %rax,%rdi
    ae08:	e8 f3 f2 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    ae0d:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    ae11:	48 89 c7             	mov    %rax,%rdi
    ae14:	e8 e7 f2 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    ae19:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    ae20:	48 89 c7             	mov    %rax,%rdi
    ae23:	e8 68 fb ff ff       	call   a990 <_ZNSaIcEC1Ev@plt>
    ae28:	48 8d 95 a8 fe ff ff 	lea    -0x158(%rbp),%rdx
    ae2f:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    ae36:	48 8d 0d 6b b5 00 00 	lea    0xb56b(%rip),%rcx        # 163a8 <_ZN3MPIL11LOCK_SHAREDE+0x168>
    ae3d:	48 89 ce             	mov    %rcx,%rsi
    ae40:	48 89 c7             	mov    %rax,%rdi
    ae43:	e8 18 f8 ff ff       	call   a660 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EPKcRKS3_@plt>
    ae48:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    ae4f:	48 89 c7             	mov    %rax,%rdi
    ae52:	e8 09 f6 ff ff       	call   a460 <_ZNSaIcED1Ev@plt>
    ae57:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    ae5e:	48 89 c7             	mov    %rax,%rdi
    ae61:	e8 2a fb ff ff       	call   a990 <_ZNSaIcEC1Ev@plt>
    ae66:	48 8d 95 a8 fe ff ff 	lea    -0x158(%rbp),%rdx
    ae6d:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    ae71:	48 8d 0d 88 b5 00 00 	lea    0xb588(%rip),%rcx        # 16400 <_ZN3MPIL11LOCK_SHAREDE+0x1c0>
    ae78:	48 89 ce             	mov    %rcx,%rsi
    ae7b:	48 89 c7             	mov    %rax,%rdi
    ae7e:	e8 dd f7 ff ff       	call   a660 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EPKcRKS3_@plt>
    ae83:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    ae8a:	48 89 c7             	mov    %rax,%rdi
    ae8d:	e8 ce f5 ff ff       	call   a460 <_ZNSaIcED1Ev@plt>
    ae92:	48 c7 85 98 fe ff ff 	movq   $0x0,-0x168(%rbp)
    ae99:	00 00 00 00 
    ae9d:	c7 85 7c fe ff ff ff 	movl   $0xffffffff,-0x184(%rbp)
    aea4:	ff ff ff 
    aea7:	c7 85 80 fe ff ff ff 	movl   $0xffffffff,-0x180(%rbp)
    aeae:	ff ff ff 
    aeb1:	c7 85 84 fe ff ff ff 	movl   $0xffffffff,-0x17c(%rbp)
    aeb8:	ff ff ff 
    aebb:	c7 85 88 fe ff ff ff 	movl   $0xffffffff,-0x178(%rbp)
    aec2:	ff ff ff 
    aec5:	48 8b 85 d0 fe ff ff 	mov    -0x130(%rbp),%rax
    aecc:	48 8d 15 4d 2f 01 00 	lea    0x12f4d(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    aed3:	48 39 d0             	cmp    %rdx,%rax
    aed6:	74 19                	je     aef1 <main+0x308>
    aed8:	48 8b 85 d0 fe ff ff 	mov    -0x130(%rbp),%rax
    aedf:	48 8d 95 7c fe ff ff 	lea    -0x184(%rbp),%rdx
    aee6:	48 89 d6             	mov    %rdx,%rsi
    aee9:	48 89 c7             	mov    %rax,%rdi
    aeec:	e8 ef f5 ff ff       	call   a4e0 <MPI_Comm_size@plt>
    aef1:	48 8b 85 d8 fe ff ff 	mov    -0x128(%rbp),%rax
    aef8:	48 8d 15 21 2f 01 00 	lea    0x12f21(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    aeff:	48 39 d0             	cmp    %rdx,%rax
    af02:	74 19                	je     af1d <main+0x334>
    af04:	48 8b 85 d8 fe ff ff 	mov    -0x128(%rbp),%rax
    af0b:	48 8d 95 80 fe ff ff 	lea    -0x180(%rbp),%rdx
    af12:	48 89 d6             	mov    %rdx,%rsi
    af15:	48 89 c7             	mov    %rax,%rdi
    af18:	e8 c3 f5 ff ff       	call   a4e0 <MPI_Comm_size@plt>
    af1d:	48 8b 85 e0 fe ff ff 	mov    -0x120(%rbp),%rax
    af24:	48 8d 15 f5 2e 01 00 	lea    0x12ef5(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    af2b:	48 39 d0             	cmp    %rdx,%rax
    af2e:	74 19                	je     af49 <main+0x360>
    af30:	48 8b 85 e0 fe ff ff 	mov    -0x120(%rbp),%rax
    af37:	48 8d 95 84 fe ff ff 	lea    -0x17c(%rbp),%rdx
    af3e:	48 89 d6             	mov    %rdx,%rsi
    af41:	48 89 c7             	mov    %rax,%rdi
    af44:	e8 97 f5 ff ff       	call   a4e0 <MPI_Comm_size@plt>
    af49:	48 8b 85 e8 fe ff ff 	mov    -0x118(%rbp),%rax
    af50:	48 8d 15 c9 2e 01 00 	lea    0x12ec9(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    af57:	48 39 d0             	cmp    %rdx,%rax
    af5a:	74 19                	je     af75 <main+0x38c>
    af5c:	48 8b 85 e8 fe ff ff 	mov    -0x118(%rbp),%rax
    af63:	48 8d 95 88 fe ff ff 	lea    -0x178(%rbp),%rdx
    af6a:	48 89 d6             	mov    %rdx,%rsi
    af6d:	48 89 c7             	mov    %rax,%rdi
    af70:	e8 6b f5 ff ff       	call   a4e0 <MPI_Comm_size@plt>
    af75:	8b 9d 88 fe ff ff    	mov    -0x178(%rbp),%ebx
    af7b:	4c 8b 9d e8 fe ff ff 	mov    -0x118(%rbp),%r11
    af82:	44 8b 95 84 fe ff ff 	mov    -0x17c(%rbp),%r10d
    af89:	4c 8b 8d e0 fe ff ff 	mov    -0x120(%rbp),%r9
    af90:	44 8b 85 80 fe ff ff 	mov    -0x180(%rbp),%r8d
    af97:	48 8b bd d8 fe ff ff 	mov    -0x128(%rbp),%rdi
    af9e:	8b b5 7c fe ff ff    	mov    -0x184(%rbp),%esi
    afa4:	4c 8b ad d0 fe ff ff 	mov    -0x130(%rbp),%r13
    afab:	44 8b a5 c4 fe ff ff 	mov    -0x13c(%rbp),%r12d
    afb2:	8b 8d c0 fe ff ff    	mov    -0x140(%rbp),%ecx
    afb8:	8b 95 cc fe ff ff    	mov    -0x134(%rbp),%edx
    afbe:	8b 85 c8 fe ff ff    	mov    -0x138(%rbp),%eax
    afc4:	48 83 ec 08          	sub    $0x8,%rsp
    afc8:	53                   	push   %rbx
    afc9:	41 53                	push   %r11
    afcb:	41 52                	push   %r10
    afcd:	41 51                	push   %r9
    afcf:	41 50                	push   %r8
    afd1:	57                   	push   %rdi
    afd2:	56                   	push   %rsi
    afd3:	4d 89 e9             	mov    %r13,%r9
    afd6:	45 89 e0             	mov    %r12d,%r8d
    afd9:	89 c6                	mov    %eax,%esi
    afdb:	48 8d 05 76 b4 00 00 	lea    0xb476(%rip),%rax        # 16458 <_ZN3MPIL11LOCK_SHAREDE+0x218>
    afe2:	48 89 c7             	mov    %rax,%rdi
    afe5:	b8 00 00 00 00       	mov    $0x0,%eax
    afea:	e8 d1 ed ff ff       	call   9dc0 <printf@plt>
    afef:	48 83 c4 40          	add    $0x40,%rsp
    aff3:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    affa:	48 85 c0             	test   %rax,%rax
    affd:	74 23                	je     b022 <main+0x439>
    afff:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b006:	48 8b 00             	mov    (%rax),%rax
    b009:	48 89 c6             	mov    %rax,%rsi
    b00c:	48 8d 05 79 b4 00 00 	lea    0xb479(%rip),%rax        # 1648c <_ZN3MPIL11LOCK_SHAREDE+0x24c>
    b013:	48 89 c7             	mov    %rax,%rdi
    b016:	b8 00 00 00 00       	mov    $0x0,%eax
    b01b:	e8 a0 ed ff ff       	call   9dc0 <printf@plt>
    b020:	eb 0f                	jmp    b031 <main+0x448>
    b022:	48 8d 05 7c b4 00 00 	lea    0xb47c(%rip),%rax        # 164a5 <_ZN3MPIL11LOCK_SHAREDE+0x265>
    b029:	48 89 c7             	mov    %rax,%rdi
    b02c:	e8 8f f7 ff ff       	call   a7c0 <puts@plt>
    b031:	48 8b 85 b0 fe ff ff 	mov    -0x150(%rbp),%rax
    b038:	8b 00                	mov    (%rax),%eax
    b03a:	85 c0                	test   %eax,%eax
    b03c:	0f 85 ae 00 00 00    	jne    b0f0 <main+0x507>
    b042:	48 8b 85 b8 fe ff ff 	mov    -0x148(%rbp),%rax
    b049:	8b 00                	mov    (%rax),%eax
    b04b:	85 c0                	test   %eax,%eax
    b04d:	0f 85 9d 00 00 00    	jne    b0f0 <main+0x507>
    b053:	48 8d 95 f0 fe ff ff 	lea    -0x110(%rbp),%rdx
    b05a:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b05e:	48 89 d6             	mov    %rdx,%rsi
    b061:	48 89 c7             	mov    %rax,%rdi
    b064:	e8 57 ef ff ff       	call   9fc0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
    b069:	48 8d 95 98 fe ff ff 	lea    -0x168(%rbp),%rdx
    b070:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b074:	48 89 d6             	mov    %rdx,%rsi
    b077:	48 89 c7             	mov    %rax,%rdi
    b07a:	e8 9d 69 00 00       	call   11a1c <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t>
    b07f:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b083:	48 89 c7             	mov    %rax,%rdi
    b086:	e8 75 f0 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b08b:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b092:	48 8d 15 4f 72 00 00 	lea    0x724f(%rip),%rdx        # 122e8 <_Z25gravityUpdateAccelerationP10particle_tS0_>
    b099:	48 89 d6             	mov    %rdx,%rsi
    b09c:	48 89 c7             	mov    %rax,%rdi
    b09f:	e8 f6 5b 00 00       	call   10c9a <_Z29registerUpdateAccelerationFnsP17chunk_particles_tPFvP10particle_tS2_E>
    b0a4:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b0ab:	be 09 00 00 00       	mov    $0x9,%esi
    b0b0:	48 89 c7             	mov    %rax,%rdi
    b0b3:	e8 da 7a 00 00       	call   12b92 <_Z7paddingP17chunk_particles_ti>
    b0b8:	48 8b 9d 98 fe ff ff 	mov    -0x168(%rbp),%rbx
    b0bf:	48 8d 95 10 ff ff ff 	lea    -0xf0(%rbp),%rdx
    b0c6:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b0ca:	48 89 d6             	mov    %rdx,%rsi
    b0cd:	48 89 c7             	mov    %rax,%rdi
    b0d0:	e8 eb ee ff ff       	call   9fc0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
    b0d5:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b0d9:	48 89 de             	mov    %rbx,%rsi
    b0dc:	48 89 c7             	mov    %rax,%rdi
    b0df:	e8 40 6a 00 00       	call   11b24 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t>
    b0e4:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b0e8:	48 89 c7             	mov    %rax,%rdi
    b0eb:	e8 10 f0 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b0f0:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b0f7:	48 85 c0             	test   %rax,%rax
    b0fa:	74 23                	je     b11f <main+0x536>
    b0fc:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b103:	48 8b 00             	mov    (%rax),%rax
    b106:	48 89 c6             	mov    %rax,%rsi
    b109:	48 8d 05 ad b3 00 00 	lea    0xb3ad(%rip),%rax        # 164bd <_ZN3MPIL11LOCK_SHAREDE+0x27d>
    b110:	48 89 c7             	mov    %rax,%rdi
    b113:	b8 00 00 00 00       	mov    $0x0,%eax
    b118:	e8 a3 ec ff ff       	call   9dc0 <printf@plt>
    b11d:	eb 0f                	jmp    b12e <main+0x545>
    b11f:	48 8d 05 b0 b3 00 00 	lea    0xb3b0(%rip),%rax        # 164d6 <_ZN3MPIL11LOCK_SHAREDE+0x296>
    b126:	48 89 c7             	mov    %rax,%rdi
    b129:	e8 92 f6 ff ff       	call   a7c0 <puts@plt>
    b12e:	48 c7 85 a0 fe ff ff 	movq   $0x0,-0x160(%rbp)
    b135:	00 00 00 00 
    b139:	48 c7 85 a8 fe ff ff 	movq   $0x0,-0x158(%rbp)
    b140:	00 00 00 00 
    b144:	48 8b 95 98 fe ff ff 	mov    -0x168(%rbp),%rdx
    b14b:	48 8d 8d c0 fe ff ff 	lea    -0x140(%rbp),%rcx
    b152:	48 8d b5 a8 fe ff ff 	lea    -0x158(%rbp),%rsi
    b159:	48 8d 85 a0 fe ff ff 	lea    -0x160(%rbp),%rax
    b160:	48 89 c7             	mov    %rax,%rdi
    b163:	e8 ef 8d 00 00       	call   13f57 <_Z7scatterPP17chunk_particles_tS1_S0_P10topology_t>
    b168:	c7 85 8c fe ff ff 00 	movl   $0x0,-0x174(%rbp)
    b16f:	00 00 00 
    b172:	eb 33                	jmp    b1a7 <main+0x5be>
    b174:	48 8b 8d a8 fe ff ff 	mov    -0x158(%rbp),%rcx
    b17b:	48 8b 85 a0 fe ff ff 	mov    -0x160(%rbp),%rax
    b182:	48 8b 35 b7 b3 00 00 	mov    0xb3b7(%rip),%rsi        # 16540 <_ZN3MPIL11LOCK_SHAREDE+0x300>
    b189:	48 8d 95 c0 fe ff ff 	lea    -0x140(%rbp),%rdx
    b190:	66 48 0f 6e c6       	movq   %rsi,%xmm0
    b195:	48 89 ce             	mov    %rcx,%rsi
    b198:	48 89 c7             	mov    %rax,%rdi
    b19b:	e8 f9 9a 00 00       	call   14c99 <_Z6updateP17chunk_particles_tS0_P10topology_td>
    b1a0:	83 85 8c fe ff ff 01 	addl   $0x1,-0x174(%rbp)
    b1a7:	83 bd 8c fe ff ff 63 	cmpl   $0x63,-0x174(%rbp)
    b1ae:	7e c4                	jle    b174 <main+0x58b>
    b1b0:	48 8b 9d a0 fe ff ff 	mov    -0x160(%rbp),%rbx
    b1b7:	48 8d 95 30 ff ff ff 	lea    -0xd0(%rbp),%rdx
    b1be:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b1c2:	48 89 d6             	mov    %rdx,%rsi
    b1c5:	48 89 c7             	mov    %rax,%rdi
    b1c8:	e8 f3 ed ff ff       	call   9fc0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
    b1cd:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b1d1:	48 89 de             	mov    %rbx,%rsi
    b1d4:	48 89 c7             	mov    %rax,%rdi
    b1d7:	e8 48 69 00 00       	call   11b24 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t>
    b1dc:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b1e0:	48 89 c7             	mov    %rax,%rdi
    b1e3:	e8 18 ef ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b1e8:	48 8b 9d a8 fe ff ff 	mov    -0x158(%rbp),%rbx
    b1ef:	48 8d 95 50 ff ff ff 	lea    -0xb0(%rbp),%rdx
    b1f6:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b1fa:	48 89 d6             	mov    %rdx,%rsi
    b1fd:	48 89 c7             	mov    %rax,%rdi
    b200:	e8 bb ed ff ff       	call   9fc0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
    b205:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b209:	48 89 de             	mov    %rbx,%rsi
    b20c:	48 89 c7             	mov    %rax,%rdi
    b20f:	e8 10 69 00 00       	call   11b24 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t>
    b214:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b218:	48 89 c7             	mov    %rax,%rdi
    b21b:	e8 e0 ee ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b220:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b227:	48 85 c0             	test   %rax,%rax
    b22a:	74 23                	je     b24f <main+0x666>
    b22c:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b233:	48 8b 00             	mov    (%rax),%rax
    b236:	48 89 c6             	mov    %rax,%rsi
    b239:	48 8d 05 ae b2 00 00 	lea    0xb2ae(%rip),%rax        # 164ee <_ZN3MPIL11LOCK_SHAREDE+0x2ae>
    b240:	48 89 c7             	mov    %rax,%rdi
    b243:	b8 00 00 00 00       	mov    $0x0,%eax
    b248:	e8 73 eb ff ff       	call   9dc0 <printf@plt>
    b24d:	eb 14                	jmp    b263 <main+0x67a>
    b24f:	48 8d 05 ab b2 00 00 	lea    0xb2ab(%rip),%rax        # 16501 <_ZN3MPIL11LOCK_SHAREDE+0x2c1>
    b256:	48 89 c7             	mov    %rax,%rdi
    b259:	b8 00 00 00 00       	mov    $0x0,%eax
    b25e:	e8 5d eb ff ff       	call   9dc0 <printf@plt>
    b263:	48 8b 95 98 fe ff ff 	mov    -0x168(%rbp),%rdx
    b26a:	48 8b b5 a8 fe ff ff 	mov    -0x158(%rbp),%rsi
    b271:	48 8b 85 a0 fe ff ff 	mov    -0x160(%rbp),%rax
    b278:	48 8d 8d c0 fe ff ff 	lea    -0x140(%rbp),%rcx
    b27f:	48 89 c7             	mov    %rax,%rdi
    b282:	e8 c1 95 00 00       	call   14848 <_Z6gatherP17chunk_particles_tS0_S0_P10topology_t>
    b287:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b28e:	48 85 c0             	test   %rax,%rax
    b291:	74 23                	je     b2b6 <main+0x6cd>
    b293:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b29a:	48 8b 00             	mov    (%rax),%rax
    b29d:	48 89 c6             	mov    %rax,%rsi
    b2a0:	48 8d 05 6d b2 00 00 	lea    0xb26d(%rip),%rax        # 16514 <_ZN3MPIL11LOCK_SHAREDE+0x2d4>
    b2a7:	48 89 c7             	mov    %rax,%rdi
    b2aa:	b8 00 00 00 00       	mov    $0x0,%eax
    b2af:	e8 0c eb ff ff       	call   9dc0 <printf@plt>
    b2b4:	eb 0f                	jmp    b2c5 <main+0x6dc>
    b2b6:	48 8d 05 6b b2 00 00 	lea    0xb26b(%rip),%rax        # 16528 <_ZN3MPIL11LOCK_SHAREDE+0x2e8>
    b2bd:	48 89 c7             	mov    %rax,%rdi
    b2c0:	e8 fb f4 ff ff       	call   a7c0 <puts@plt>
    b2c5:	48 8b 85 b0 fe ff ff 	mov    -0x150(%rbp),%rax
    b2cc:	8b 00                	mov    (%rax),%eax
    b2ce:	85 c0                	test   %eax,%eax
    b2d0:	0f 85 89 00 00 00    	jne    b35f <main+0x776>
    b2d6:	48 8b 85 b8 fe ff ff 	mov    -0x148(%rbp),%rax
    b2dd:	8b 00                	mov    (%rax),%eax
    b2df:	85 c0                	test   %eax,%eax
    b2e1:	75 7c                	jne    b35f <main+0x776>
    b2e3:	48 8b 9d 98 fe ff ff 	mov    -0x168(%rbp),%rbx
    b2ea:	48 8d 95 70 ff ff ff 	lea    -0x90(%rbp),%rdx
    b2f1:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b2f5:	48 89 d6             	mov    %rdx,%rsi
    b2f8:	48 89 c7             	mov    %rax,%rdi
    b2fb:	e8 c0 ec ff ff       	call   9fc0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
    b300:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b304:	48 89 de             	mov    %rbx,%rsi
    b307:	48 89 c7             	mov    %rax,%rdi
    b30a:	e8 15 68 00 00       	call   11b24 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t>
    b30f:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b313:	48 89 c7             	mov    %rax,%rdi
    b316:	e8 e5 ed ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b31b:	48 8b 85 98 fe ff ff 	mov    -0x168(%rbp),%rax
    b322:	48 89 c7             	mov    %rax,%rdi
    b325:	e8 8f 7c 00 00       	call   12fb9 <_Z6shrinkP17chunk_particles_t>
    b32a:	48 8b 9d 98 fe ff ff 	mov    -0x168(%rbp),%rbx
    b331:	48 8d 55 90          	lea    -0x70(%rbp),%rdx
    b335:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b339:	48 89 d6             	mov    %rdx,%rsi
    b33c:	48 89 c7             	mov    %rax,%rdi
    b33f:	e8 7c ec ff ff       	call   9fc0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
    b344:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b348:	48 89 de             	mov    %rbx,%rsi
    b34b:	48 89 c7             	mov    %rax,%rdi
    b34e:	e8 d1 67 00 00       	call   11b24 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t>
    b353:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b357:	48 89 c7             	mov    %rax,%rdi
    b35a:	e8 a1 ed ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b35f:	e8 ec ee ff ff       	call   a250 <MPI_Finalize@plt>
    b364:	bb 00 00 00 00       	mov    $0x0,%ebx
    b369:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    b36d:	48 89 c7             	mov    %rax,%rdi
    b370:	e8 8b ed ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b375:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    b37c:	48 89 c7             	mov    %rax,%rdi
    b37f:	e8 7c ed ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b384:	48 8d 85 50 ff ff ff 	lea    -0xb0(%rbp),%rax
    b38b:	48 89 c7             	mov    %rax,%rdi
    b38e:	e8 6d ed ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b393:	48 8d 85 30 ff ff ff 	lea    -0xd0(%rbp),%rax
    b39a:	48 89 c7             	mov    %rax,%rdi
    b39d:	e8 5e ed ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b3a2:	48 8d 85 10 ff ff ff 	lea    -0xf0(%rbp),%rax
    b3a9:	48 89 c7             	mov    %rax,%rdi
    b3ac:	e8 4f ed ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b3b1:	48 8d 85 f0 fe ff ff 	lea    -0x110(%rbp),%rax
    b3b8:	48 89 c7             	mov    %rax,%rdi
    b3bb:	e8 40 ed ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b3c0:	89 d8                	mov    %ebx,%eax
    b3c2:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
    b3c6:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    b3cd:	00 00 
    b3cf:	0f 84 ff 01 00 00    	je     b5d4 <main+0x9eb>
    b3d5:	e9 f5 01 00 00       	jmp    b5cf <main+0x9e6>
    b3da:	f3 0f 1e fa          	endbr64 
    b3de:	48 89 c3             	mov    %rax,%rbx
    b3e1:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    b3e8:	48 89 c7             	mov    %rax,%rdi
    b3eb:	e8 70 f0 ff ff       	call   a460 <_ZNSaIcED1Ev@plt>
    b3f0:	48 89 d8             	mov    %rbx,%rax
    b3f3:	48 89 c7             	mov    %rax,%rdi
    b3f6:	e8 85 f5 ff ff       	call   a980 <_Unwind_Resume@plt>
    b3fb:	f3 0f 1e fa          	endbr64 
    b3ff:	48 89 c3             	mov    %rax,%rbx
    b402:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    b409:	48 89 c7             	mov    %rax,%rdi
    b40c:	e8 4f f0 ff ff       	call   a460 <_ZNSaIcED1Ev@plt>
    b411:	e9 9f 01 00 00       	jmp    b5b5 <main+0x9cc>
    b416:	f3 0f 1e fa          	endbr64 
    b41a:	48 89 c3             	mov    %rax,%rbx
    b41d:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    b421:	48 89 c7             	mov    %rax,%rdi
    b424:	e8 d7 ec ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b429:	eb 07                	jmp    b432 <main+0x849>
    b42b:	f3 0f 1e fa          	endbr64 
    b42f:	48 89 c3             	mov    %rax,%rbx
    b432:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    b439:	48 89 c7             	mov    %rax,%rdi
    b43c:	e8 bf ec ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b441:	eb 07                	jmp    b44a <main+0x861>
    b443:	f3 0f 1e fa          	endbr64 
    b447:	48 89 c3             	mov    %rax,%rbx
    b44a:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b44e:	48 89 c7             	mov    %rax,%rdi
    b451:	e8 aa ec ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b456:	e9 4b 01 00 00       	jmp    b5a6 <main+0x9bd>
    b45b:	f3 0f 1e fa          	endbr64 
    b45f:	48 89 c3             	mov    %rax,%rbx
    b462:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    b466:	48 89 c7             	mov    %rax,%rdi
    b469:	e8 92 ec ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b46e:	eb 07                	jmp    b477 <main+0x88e>
    b470:	f3 0f 1e fa          	endbr64 
    b474:	48 89 c3             	mov    %rax,%rbx
    b477:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    b47e:	48 89 c7             	mov    %rax,%rdi
    b481:	e8 7a ec ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b486:	eb 07                	jmp    b48f <main+0x8a6>
    b488:	f3 0f 1e fa          	endbr64 
    b48c:	48 89 c3             	mov    %rax,%rbx
    b48f:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b493:	48 89 c7             	mov    %rax,%rdi
    b496:	e8 65 ec ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b49b:	e9 ee 00 00 00       	jmp    b58e <main+0x9a5>
    b4a0:	f3 0f 1e fa          	endbr64 
    b4a4:	48 89 c3             	mov    %rax,%rbx
    b4a7:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    b4ae:	48 89 c7             	mov    %rax,%rdi
    b4b1:	e8 aa ef ff ff       	call   a460 <_ZNSaIcED1Ev@plt>
    b4b6:	e9 bb 00 00 00       	jmp    b576 <main+0x98d>
    b4bb:	f3 0f 1e fa          	endbr64 
    b4bf:	48 89 c3             	mov    %rax,%rbx
    b4c2:	48 8d 85 a8 fe ff ff 	lea    -0x158(%rbp),%rax
    b4c9:	48 89 c7             	mov    %rax,%rdi
    b4cc:	e8 8f ef ff ff       	call   a460 <_ZNSaIcED1Ev@plt>
    b4d1:	e9 91 00 00 00       	jmp    b567 <main+0x97e>
    b4d6:	f3 0f 1e fa          	endbr64 
    b4da:	48 89 c3             	mov    %rax,%rbx
    b4dd:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b4e1:	48 89 c7             	mov    %rax,%rdi
    b4e4:	e8 17 ec ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b4e9:	eb 70                	jmp    b55b <main+0x972>
    b4eb:	f3 0f 1e fa          	endbr64 
    b4ef:	48 89 c3             	mov    %rax,%rbx
    b4f2:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b4f6:	48 89 c7             	mov    %rax,%rdi
    b4f9:	e8 02 ec ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b4fe:	eb 5b                	jmp    b55b <main+0x972>
    b500:	f3 0f 1e fa          	endbr64 
    b504:	48 89 c3             	mov    %rax,%rbx
    b507:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b50b:	48 89 c7             	mov    %rax,%rdi
    b50e:	e8 ed eb ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b513:	eb 46                	jmp    b55b <main+0x972>
    b515:	f3 0f 1e fa          	endbr64 
    b519:	48 89 c3             	mov    %rax,%rbx
    b51c:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b520:	48 89 c7             	mov    %rax,%rdi
    b523:	e8 d8 eb ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b528:	eb 31                	jmp    b55b <main+0x972>
    b52a:	f3 0f 1e fa          	endbr64 
    b52e:	48 89 c3             	mov    %rax,%rbx
    b531:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b535:	48 89 c7             	mov    %rax,%rdi
    b538:	e8 c3 eb ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b53d:	eb 1c                	jmp    b55b <main+0x972>
    b53f:	f3 0f 1e fa          	endbr64 
    b543:	48 89 c3             	mov    %rax,%rbx
    b546:	48 8d 45 b0          	lea    -0x50(%rbp),%rax
    b54a:	48 89 c7             	mov    %rax,%rdi
    b54d:	e8 ae eb ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b552:	eb 07                	jmp    b55b <main+0x972>
    b554:	f3 0f 1e fa          	endbr64 
    b558:	48 89 c3             	mov    %rax,%rbx
    b55b:	48 8d 45 90          	lea    -0x70(%rbp),%rax
    b55f:	48 89 c7             	mov    %rax,%rdi
    b562:	e8 99 eb ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b567:	48 8d 85 70 ff ff ff 	lea    -0x90(%rbp),%rax
    b56e:	48 89 c7             	mov    %rax,%rdi
    b571:	e8 8a eb ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b576:	48 8d 85 50 ff ff ff 	lea    -0xb0(%rbp),%rax
    b57d:	48 89 c7             	mov    %rax,%rdi
    b580:	e8 7b eb ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b585:	eb 07                	jmp    b58e <main+0x9a5>
    b587:	f3 0f 1e fa          	endbr64 
    b58b:	48 89 c3             	mov    %rax,%rbx
    b58e:	48 8d 85 30 ff ff ff 	lea    -0xd0(%rbp),%rax
    b595:	48 89 c7             	mov    %rax,%rdi
    b598:	e8 63 eb ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b59d:	eb 07                	jmp    b5a6 <main+0x9bd>
    b59f:	f3 0f 1e fa          	endbr64 
    b5a3:	48 89 c3             	mov    %rax,%rbx
    b5a6:	48 8d 85 10 ff ff ff 	lea    -0xf0(%rbp),%rax
    b5ad:	48 89 c7             	mov    %rax,%rdi
    b5b0:	e8 4b eb ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b5b5:	48 8d 85 f0 fe ff ff 	lea    -0x110(%rbp),%rax
    b5bc:	48 89 c7             	mov    %rax,%rdi
    b5bf:	e8 3c eb ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b5c4:	48 89 d8             	mov    %rbx,%rax
    b5c7:	48 89 c7             	mov    %rax,%rdi
    b5ca:	e8 b1 f3 ff ff       	call   a980 <_Unwind_Resume@plt>
    b5cf:	e8 fc ee ff ff       	call   a4d0 <__stack_chk_fail@plt>
    b5d4:	48 8d 65 e8          	lea    -0x18(%rbp),%rsp
    b5d8:	5b                   	pop    %rbx
    b5d9:	41 5c                	pop    %r12
    b5db:	41 5d                	pop    %r13
    b5dd:	5d                   	pop    %rbp
    b5de:	c3                   	ret    

000000000000b5df <_Z41__static_initialization_and_destruction_0ii>:
    b5df:	f3 0f 1e fa          	endbr64 
    b5e3:	55                   	push   %rbp
    b5e4:	48 89 e5             	mov    %rsp,%rbp
    b5e7:	48 83 ec 10          	sub    $0x10,%rsp
    b5eb:	89 7d fc             	mov    %edi,-0x4(%rbp)
    b5ee:	89 75 f8             	mov    %esi,-0x8(%rbp)
    b5f1:	83 7d fc 01          	cmpl   $0x1,-0x4(%rbp)
    b5f5:	75 3b                	jne    b632 <_Z41__static_initialization_and_destruction_0ii+0x53>
    b5f7:	81 7d f8 ff ff 00 00 	cmpl   $0xffff,-0x8(%rbp)
    b5fe:	75 32                	jne    b632 <_Z41__static_initialization_and_destruction_0ii+0x53>
    b600:	48 8d 05 22 2c 01 00 	lea    0x12c22(%rip),%rax        # 1e229 <_ZStL8__ioinit>
    b607:	48 89 c7             	mov    %rax,%rdi
    b60a:	e8 a1 f1 ff ff       	call   a7b0 <_ZNSt8ios_base4InitC1Ev@plt>
    b60f:	48 8d 05 f2 19 01 00 	lea    0x119f2(%rip),%rax        # 1d008 <__dso_handle>
    b616:	48 89 c2             	mov    %rax,%rdx
    b619:	48 8d 05 09 2c 01 00 	lea    0x12c09(%rip),%rax        # 1e229 <_ZStL8__ioinit>
    b620:	48 89 c6             	mov    %rax,%rsi
    b623:	48 8b 05 ce 19 01 00 	mov    0x119ce(%rip),%rax        # 1cff8 <_ZNSt8ios_base4InitD1Ev@GLIBCXX_3.4>
    b62a:	48 89 c7             	mov    %rax,%rdi
    b62d:	e8 ee ec ff ff       	call   a320 <__cxa_atexit@plt>
    b632:	90                   	nop
    b633:	c9                   	leave  
    b634:	c3                   	ret    

000000000000b635 <_GLOBAL__sub_I_main>:
    b635:	f3 0f 1e fa          	endbr64 
    b639:	55                   	push   %rbp
    b63a:	48 89 e5             	mov    %rsp,%rbp
    b63d:	be ff ff 00 00       	mov    $0xffff,%esi
    b642:	bf 01 00 00 00       	mov    $0x1,%edi
    b647:	e8 93 ff ff ff       	call   b5df <_Z41__static_initialization_and_destruction_0ii>
    b64c:	5d                   	pop    %rbp
    b64d:	c3                   	ret    

000000000000b64e <_ZNSt7__cxx119to_stringEi>:
    b64e:	f3 0f 1e fa          	endbr64 
    b652:	55                   	push   %rbp
    b653:	48 89 e5             	mov    %rsp,%rbp
    b656:	53                   	push   %rbx
    b657:	48 83 ec 38          	sub    $0x38,%rsp
    b65b:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
    b65f:	89 75 c4             	mov    %esi,-0x3c(%rbp)
    b662:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    b669:	00 00 
    b66b:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    b66f:	31 c0                	xor    %eax,%eax
    b671:	8b 45 c4             	mov    -0x3c(%rbp),%eax
    b674:	c1 e8 1f             	shr    $0x1f,%eax
    b677:	88 45 df             	mov    %al,-0x21(%rbp)
    b67a:	80 7d df 00          	cmpb   $0x0,-0x21(%rbp)
    b67e:	74 07                	je     b687 <_ZNSt7__cxx119to_stringEi+0x39>
    b680:	8b 45 c4             	mov    -0x3c(%rbp),%eax
    b683:	f7 d8                	neg    %eax
    b685:	eb 03                	jmp    b68a <_ZNSt7__cxx119to_stringEi+0x3c>
    b687:	8b 45 c4             	mov    -0x3c(%rbp),%eax
    b68a:	89 45 e0             	mov    %eax,-0x20(%rbp)
    b68d:	8b 45 e0             	mov    -0x20(%rbp),%eax
    b690:	be 0a 00 00 00       	mov    $0xa,%esi
    b695:	89 c7                	mov    %eax,%edi
    b697:	e8 c1 00 00 00       	call   b75d <_ZNSt8__detail14__to_chars_lenIjEEjT_i>
    b69c:	89 45 e4             	mov    %eax,-0x1c(%rbp)
    b69f:	48 8d 45 de          	lea    -0x22(%rbp),%rax
    b6a3:	48 89 c7             	mov    %rax,%rdi
    b6a6:	e8 e5 f2 ff ff       	call   a990 <_ZNSaIcEC1Ev@plt>
    b6ab:	0f b6 55 df          	movzbl -0x21(%rbp),%edx
    b6af:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    b6b2:	01 d0                	add    %edx,%eax
    b6b4:	89 c6                	mov    %eax,%esi
    b6b6:	48 8d 55 de          	lea    -0x22(%rbp),%rdx
    b6ba:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    b6be:	48 89 d1             	mov    %rdx,%rcx
    b6c1:	ba 2d 00 00 00       	mov    $0x2d,%edx
    b6c6:	48 89 c7             	mov    %rax,%rdi
    b6c9:	e8 72 ea ff ff       	call   a140 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EmcRKS3_@plt>
    b6ce:	48 8d 45 de          	lea    -0x22(%rbp),%rax
    b6d2:	48 89 c7             	mov    %rax,%rdi
    b6d5:	e8 86 ed ff ff       	call   a460 <_ZNSaIcED1Ev@plt>
    b6da:	0f b6 55 df          	movzbl -0x21(%rbp),%edx
    b6de:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    b6e2:	48 89 d6             	mov    %rdx,%rsi
    b6e5:	48 89 c7             	mov    %rax,%rdi
    b6e8:	e8 03 f4 ff ff       	call   aaf0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEixEm@plt>
    b6ed:	48 89 c1             	mov    %rax,%rcx
    b6f0:	8b 55 e0             	mov    -0x20(%rbp),%edx
    b6f3:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    b6f6:	89 c6                	mov    %eax,%esi
    b6f8:	48 89 cf             	mov    %rcx,%rdi
    b6fb:	e8 1e 4f 00 00       	call   1061e <_ZNSt8__detail18__to_chars_10_implIjEEvPcjT_>
    b700:	90                   	nop
    b701:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    b705:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    b70c:	00 00 
    b70e:	74 43                	je     b753 <_ZNSt7__cxx119to_stringEi+0x105>
    b710:	eb 3c                	jmp    b74e <_ZNSt7__cxx119to_stringEi+0x100>
    b712:	f3 0f 1e fa          	endbr64 
    b716:	48 89 c3             	mov    %rax,%rbx
    b719:	48 8d 45 de          	lea    -0x22(%rbp),%rax
    b71d:	48 89 c7             	mov    %rax,%rdi
    b720:	e8 3b ed ff ff       	call   a460 <_ZNSaIcED1Ev@plt>
    b725:	48 89 d8             	mov    %rbx,%rax
    b728:	48 89 c7             	mov    %rax,%rdi
    b72b:	e8 50 f2 ff ff       	call   a980 <_Unwind_Resume@plt>
    b730:	f3 0f 1e fa          	endbr64 
    b734:	48 89 c3             	mov    %rax,%rbx
    b737:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    b73b:	48 89 c7             	mov    %rax,%rdi
    b73e:	e8 bd e9 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
    b743:	48 89 d8             	mov    %rbx,%rax
    b746:	48 89 c7             	mov    %rax,%rdi
    b749:	e8 32 f2 ff ff       	call   a980 <_Unwind_Resume@plt>
    b74e:	e8 7d ed ff ff       	call   a4d0 <__stack_chk_fail@plt>
    b753:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    b757:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    b75b:	c9                   	leave  
    b75c:	c3                   	ret    

000000000000b75d <_ZNSt8__detail14__to_chars_lenIjEEjT_i>:
    b75d:	f3 0f 1e fa          	endbr64 
    b761:	55                   	push   %rbp
    b762:	48 89 e5             	mov    %rsp,%rbp
    b765:	89 7d dc             	mov    %edi,-0x24(%rbp)
    b768:	89 75 d8             	mov    %esi,-0x28(%rbp)
    b76b:	c7 45 ec 01 00 00 00 	movl   $0x1,-0x14(%rbp)
    b772:	8b 45 d8             	mov    -0x28(%rbp),%eax
    b775:	0f af c0             	imul   %eax,%eax
    b778:	89 45 f0             	mov    %eax,-0x10(%rbp)
    b77b:	8b 45 d8             	mov    -0x28(%rbp),%eax
    b77e:	8b 55 f0             	mov    -0x10(%rbp),%edx
    b781:	0f af c2             	imul   %edx,%eax
    b784:	89 45 f4             	mov    %eax,-0xc(%rbp)
    b787:	8b 45 d8             	mov    -0x28(%rbp),%eax
    b78a:	0f af 45 f4          	imul   -0xc(%rbp),%eax
    b78e:	89 c0                	mov    %eax,%eax
    b790:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    b794:	8b 45 d8             	mov    -0x28(%rbp),%eax
    b797:	39 45 dc             	cmp    %eax,-0x24(%rbp)
    b79a:	73 05                	jae    b7a1 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x44>
    b79c:	8b 45 ec             	mov    -0x14(%rbp),%eax
    b79f:	eb 46                	jmp    b7e7 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x8a>
    b7a1:	8b 45 dc             	mov    -0x24(%rbp),%eax
    b7a4:	3b 45 f0             	cmp    -0x10(%rbp),%eax
    b7a7:	73 08                	jae    b7b1 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x54>
    b7a9:	8b 45 ec             	mov    -0x14(%rbp),%eax
    b7ac:	83 c0 01             	add    $0x1,%eax
    b7af:	eb 36                	jmp    b7e7 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x8a>
    b7b1:	8b 45 dc             	mov    -0x24(%rbp),%eax
    b7b4:	3b 45 f4             	cmp    -0xc(%rbp),%eax
    b7b7:	73 08                	jae    b7c1 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x64>
    b7b9:	8b 45 ec             	mov    -0x14(%rbp),%eax
    b7bc:	83 c0 02             	add    $0x2,%eax
    b7bf:	eb 26                	jmp    b7e7 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x8a>
    b7c1:	8b 45 dc             	mov    -0x24(%rbp),%eax
    b7c4:	48 39 45 f8          	cmp    %rax,-0x8(%rbp)
    b7c8:	76 08                	jbe    b7d2 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x75>
    b7ca:	8b 45 ec             	mov    -0x14(%rbp),%eax
    b7cd:	83 c0 03             	add    $0x3,%eax
    b7d0:	eb 15                	jmp    b7e7 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x8a>
    b7d2:	8b 45 dc             	mov    -0x24(%rbp),%eax
    b7d5:	ba 00 00 00 00       	mov    $0x0,%edx
    b7da:	48 f7 75 f8          	divq   -0x8(%rbp)
    b7de:	89 45 dc             	mov    %eax,-0x24(%rbp)
    b7e1:	83 45 ec 04          	addl   $0x4,-0x14(%rbp)
    b7e5:	eb ad                	jmp    b794 <_ZNSt8__detail14__to_chars_lenIjEEjT_i+0x37>
    b7e7:	5d                   	pop    %rbp
    b7e8:	c3                   	ret    
    b7e9:	90                   	nop

000000000000b7ea <_ZNKSt17integral_constantIbLb1EEcvbEv>:
    b7ea:	f3 0f 1e fa          	endbr64 
    b7ee:	55                   	push   %rbp
    b7ef:	48 89 e5             	mov    %rsp,%rbp
    b7f2:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b7f6:	b8 01 00 00 00       	mov    $0x1,%eax
    b7fb:	5d                   	pop    %rbp
    b7fc:	c3                   	ret    
    b7fd:	90                   	nop

000000000000b7fe <_ZN3MPI8DatatypeD1Ev>:
    b7fe:	f3 0f 1e fa          	endbr64 
    b802:	55                   	push   %rbp
    b803:	48 89 e5             	mov    %rsp,%rbp
    b806:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b80a:	48 8d 15 6f 0c 01 00 	lea    0x10c6f(%rip),%rdx        # 1c480 <_ZTVN3MPI8DatatypeE+0x10>
    b811:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b815:	48 89 10             	mov    %rdx,(%rax)
    b818:	90                   	nop
    b819:	5d                   	pop    %rbp
    b81a:	c3                   	ret    
    b81b:	90                   	nop

000000000000b81c <_ZN3MPI8DatatypeD0Ev>:
    b81c:	f3 0f 1e fa          	endbr64 
    b820:	55                   	push   %rbp
    b821:	48 89 e5             	mov    %rsp,%rbp
    b824:	48 83 ec 10          	sub    $0x10,%rsp
    b828:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b82c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b830:	48 89 c7             	mov    %rax,%rdi
    b833:	e8 c6 ff ff ff       	call   b7fe <_ZN3MPI8DatatypeD1Ev>
    b838:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b83c:	be 10 00 00 00       	mov    $0x10,%esi
    b841:	48 89 c7             	mov    %rax,%rdi
    b844:	e8 a7 eb ff ff       	call   a3f0 <_ZdlPvm@plt>
    b849:	c9                   	leave  
    b84a:	c3                   	ret    
    b84b:	90                   	nop

000000000000b84c <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>:
    b84c:	f3 0f 1e fa          	endbr64 
    b850:	55                   	push   %rbp
    b851:	48 89 e5             	mov    %rsp,%rbp
    b854:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b858:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    b85c:	48 8d 15 1d 0c 01 00 	lea    0x10c1d(%rip),%rdx        # 1c480 <_ZTVN3MPI8DatatypeE+0x10>
    b863:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b867:	48 89 10             	mov    %rdx,(%rax)
    b86a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b86e:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    b872:	48 89 50 08          	mov    %rdx,0x8(%rax)
    b876:	90                   	nop
    b877:	5d                   	pop    %rbp
    b878:	c3                   	ret    
    b879:	90                   	nop

000000000000b87a <_ZN3MPI8DatatypeaSERKP15ompi_datatype_t>:
    b87a:	f3 0f 1e fa          	endbr64 
    b87e:	55                   	push   %rbp
    b87f:	48 89 e5             	mov    %rsp,%rbp
    b882:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b886:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    b88a:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    b88e:	48 8b 10             	mov    (%rax),%rdx
    b891:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b895:	48 89 50 08          	mov    %rdx,0x8(%rax)
    b899:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b89d:	5d                   	pop    %rbp
    b89e:	c3                   	ret    
    b89f:	90                   	nop

000000000000b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>:
    b8a0:	f3 0f 1e fa          	endbr64 
    b8a4:	55                   	push   %rbp
    b8a5:	48 89 e5             	mov    %rsp,%rbp
    b8a8:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b8ac:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b8b0:	48 8b 40 08          	mov    0x8(%rax),%rax
    b8b4:	5d                   	pop    %rbp
    b8b5:	c3                   	ret    

000000000000b8b6 <_ZN3MPI6StatusD1Ev>:
    b8b6:	f3 0f 1e fa          	endbr64 
    b8ba:	55                   	push   %rbp
    b8bb:	48 89 e5             	mov    %rsp,%rbp
    b8be:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b8c2:	48 8d 15 9f fd 00 00 	lea    0xfd9f(%rip),%rdx        # 1b668 <_ZTVN3MPI6StatusE+0x10>
    b8c9:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b8cd:	48 89 10             	mov    %rdx,(%rax)
    b8d0:	90                   	nop
    b8d1:	5d                   	pop    %rbp
    b8d2:	c3                   	ret    
    b8d3:	90                   	nop

000000000000b8d4 <_ZN3MPI6StatusD0Ev>:
    b8d4:	f3 0f 1e fa          	endbr64 
    b8d8:	55                   	push   %rbp
    b8d9:	48 89 e5             	mov    %rsp,%rbp
    b8dc:	48 83 ec 10          	sub    $0x10,%rsp
    b8e0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b8e4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b8e8:	48 89 c7             	mov    %rax,%rdi
    b8eb:	e8 c6 ff ff ff       	call   b8b6 <_ZN3MPI6StatusD1Ev>
    b8f0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b8f4:	be 20 00 00 00       	mov    $0x20,%esi
    b8f9:	48 89 c7             	mov    %rax,%rdi
    b8fc:	e8 ef ea ff ff       	call   a3f0 <_ZdlPvm@plt>
    b901:	c9                   	leave  
    b902:	c3                   	ret    
    b903:	90                   	nop

000000000000b904 <_ZN3MPI6StatusaSERK20ompi_status_public_t>:
    b904:	f3 0f 1e fa          	endbr64 
    b908:	55                   	push   %rbp
    b909:	48 89 e5             	mov    %rsp,%rbp
    b90c:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b910:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    b914:	48 8b 4d f8          	mov    -0x8(%rbp),%rcx
    b918:	48 8b 75 f0          	mov    -0x10(%rbp),%rsi
    b91c:	48 8b 06             	mov    (%rsi),%rax
    b91f:	48 8b 56 08          	mov    0x8(%rsi),%rdx
    b923:	48 89 41 08          	mov    %rax,0x8(%rcx)
    b927:	48 89 51 10          	mov    %rdx,0x10(%rcx)
    b92b:	48 8b 46 10          	mov    0x10(%rsi),%rax
    b92f:	48 89 41 18          	mov    %rax,0x18(%rcx)
    b933:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b937:	5d                   	pop    %rbp
    b938:	c3                   	ret    
    b939:	90                   	nop

000000000000b93a <_ZN3MPI7RequestD1Ev>:
    b93a:	f3 0f 1e fa          	endbr64 
    b93e:	55                   	push   %rbp
    b93f:	48 89 e5             	mov    %rsp,%rbp
    b942:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b946:	48 8d 15 23 08 01 00 	lea    0x10823(%rip),%rdx        # 1c170 <_ZTVN3MPI7RequestE+0x10>
    b94d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b951:	48 89 10             	mov    %rdx,(%rax)
    b954:	90                   	nop
    b955:	5d                   	pop    %rbp
    b956:	c3                   	ret    
    b957:	90                   	nop

000000000000b958 <_ZN3MPI7RequestD0Ev>:
    b958:	f3 0f 1e fa          	endbr64 
    b95c:	55                   	push   %rbp
    b95d:	48 89 e5             	mov    %rsp,%rbp
    b960:	48 83 ec 10          	sub    $0x10,%rsp
    b964:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b968:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b96c:	48 89 c7             	mov    %rax,%rdi
    b96f:	e8 c6 ff ff ff       	call   b93a <_ZN3MPI7RequestD1Ev>
    b974:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b978:	be 10 00 00 00       	mov    $0x10,%esi
    b97d:	48 89 c7             	mov    %rax,%rdi
    b980:	e8 6b ea ff ff       	call   a3f0 <_ZdlPvm@plt>
    b985:	c9                   	leave  
    b986:	c3                   	ret    
    b987:	90                   	nop

000000000000b988 <_ZN3MPI7RequestC1EP14ompi_request_t>:
    b988:	f3 0f 1e fa          	endbr64 
    b98c:	55                   	push   %rbp
    b98d:	48 89 e5             	mov    %rsp,%rbp
    b990:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b994:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    b998:	48 8d 15 d1 07 01 00 	lea    0x107d1(%rip),%rdx        # 1c170 <_ZTVN3MPI7RequestE+0x10>
    b99f:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b9a3:	48 89 10             	mov    %rdx,(%rax)
    b9a6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b9aa:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    b9ae:	48 89 50 08          	mov    %rdx,0x8(%rax)
    b9b2:	90                   	nop
    b9b3:	5d                   	pop    %rbp
    b9b4:	c3                   	ret    
    b9b5:	90                   	nop

000000000000b9b6 <_ZN3MPI8PrequestC1ERKP14ompi_request_t>:
    b9b6:	f3 0f 1e fa          	endbr64 
    b9ba:	55                   	push   %rbp
    b9bb:	48 89 e5             	mov    %rsp,%rbp
    b9be:	48 83 ec 10          	sub    $0x10,%rsp
    b9c2:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    b9c6:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    b9ca:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b9ce:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    b9d2:	48 8b 12             	mov    (%rdx),%rdx
    b9d5:	48 89 d6             	mov    %rdx,%rsi
    b9d8:	48 89 c7             	mov    %rax,%rdi
    b9db:	e8 a8 ff ff ff       	call   b988 <_ZN3MPI7RequestC1EP14ompi_request_t>
    b9e0:	48 8d 15 21 07 01 00 	lea    0x10721(%rip),%rdx        # 1c108 <_ZTVN3MPI8PrequestE+0x10>
    b9e7:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    b9eb:	48 89 10             	mov    %rdx,(%rax)
    b9ee:	90                   	nop
    b9ef:	c9                   	leave  
    b9f0:	c3                   	ret    
    b9f1:	90                   	nop

000000000000b9f2 <_ZN3MPI8PrequestD1Ev>:
    b9f2:	f3 0f 1e fa          	endbr64 
    b9f6:	55                   	push   %rbp
    b9f7:	48 89 e5             	mov    %rsp,%rbp
    b9fa:	48 83 ec 10          	sub    $0x10,%rsp
    b9fe:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    ba02:	48 8d 15 ff 06 01 00 	lea    0x106ff(%rip),%rdx        # 1c108 <_ZTVN3MPI8PrequestE+0x10>
    ba09:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba0d:	48 89 10             	mov    %rdx,(%rax)
    ba10:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba14:	48 89 c7             	mov    %rax,%rdi
    ba17:	e8 1e ff ff ff       	call   b93a <_ZN3MPI7RequestD1Ev>
    ba1c:	90                   	nop
    ba1d:	c9                   	leave  
    ba1e:	c3                   	ret    
    ba1f:	90                   	nop

000000000000ba20 <_ZN3MPI8PrequestD0Ev>:
    ba20:	f3 0f 1e fa          	endbr64 
    ba24:	55                   	push   %rbp
    ba25:	48 89 e5             	mov    %rsp,%rbp
    ba28:	48 83 ec 10          	sub    $0x10,%rsp
    ba2c:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    ba30:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba34:	48 89 c7             	mov    %rax,%rdi
    ba37:	e8 b6 ff ff ff       	call   b9f2 <_ZN3MPI8PrequestD1Ev>
    ba3c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba40:	be 10 00 00 00       	mov    $0x10,%esi
    ba45:	48 89 c7             	mov    %rax,%rdi
    ba48:	e8 a3 e9 ff ff       	call   a3f0 <_ZdlPvm@plt>
    ba4d:	c9                   	leave  
    ba4e:	c3                   	ret    
    ba4f:	90                   	nop

000000000000ba50 <_ZN3MPI8GrequestD1Ev>:
    ba50:	f3 0f 1e fa          	endbr64 
    ba54:	55                   	push   %rbp
    ba55:	48 89 e5             	mov    %rsp,%rbp
    ba58:	48 83 ec 10          	sub    $0x10,%rsp
    ba5c:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    ba60:	48 8d 15 39 06 01 00 	lea    0x10639(%rip),%rdx        # 1c0a0 <_ZTVN3MPI8GrequestE+0x10>
    ba67:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba6b:	48 89 10             	mov    %rdx,(%rax)
    ba6e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba72:	48 89 c7             	mov    %rax,%rdi
    ba75:	e8 c0 fe ff ff       	call   b93a <_ZN3MPI7RequestD1Ev>
    ba7a:	90                   	nop
    ba7b:	c9                   	leave  
    ba7c:	c3                   	ret    
    ba7d:	90                   	nop

000000000000ba7e <_ZN3MPI8GrequestD0Ev>:
    ba7e:	f3 0f 1e fa          	endbr64 
    ba82:	55                   	push   %rbp
    ba83:	48 89 e5             	mov    %rsp,%rbp
    ba86:	48 83 ec 10          	sub    $0x10,%rsp
    ba8a:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    ba8e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba92:	48 89 c7             	mov    %rax,%rdi
    ba95:	e8 b6 ff ff ff       	call   ba50 <_ZN3MPI8GrequestD1Ev>
    ba9a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ba9e:	be 10 00 00 00       	mov    $0x10,%esi
    baa3:	48 89 c7             	mov    %rax,%rdi
    baa6:	e8 45 e9 ff ff       	call   a3f0 <_ZdlPvm@plt>
    baab:	c9                   	leave  
    baac:	c3                   	ret    
    baad:	90                   	nop

000000000000baae <_ZN3MPI5GroupC1EP12ompi_group_t>:
    baae:	f3 0f 1e fa          	endbr64 
    bab2:	55                   	push   %rbp
    bab3:	48 89 e5             	mov    %rsp,%rbp
    bab6:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    baba:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    babe:	48 8d 15 83 fc 00 00 	lea    0xfc83(%rip),%rdx        # 1b748 <_ZTVN3MPI5GroupE+0x10>
    bac5:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bac9:	48 89 10             	mov    %rdx,(%rax)
    bacc:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bad0:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bad4:	48 89 50 08          	mov    %rdx,0x8(%rax)
    bad8:	90                   	nop
    bad9:	5d                   	pop    %rbp
    bada:	c3                   	ret    
    badb:	90                   	nop

000000000000badc <_ZN3MPI5GroupD1Ev>:
    badc:	f3 0f 1e fa          	endbr64 
    bae0:	55                   	push   %rbp
    bae1:	48 89 e5             	mov    %rsp,%rbp
    bae4:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bae8:	48 8d 15 59 fc 00 00 	lea    0xfc59(%rip),%rdx        # 1b748 <_ZTVN3MPI5GroupE+0x10>
    baef:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    baf3:	48 89 10             	mov    %rdx,(%rax)
    baf6:	90                   	nop
    baf7:	5d                   	pop    %rbp
    baf8:	c3                   	ret    
    baf9:	90                   	nop

000000000000bafa <_ZN3MPI5GroupD0Ev>:
    bafa:	f3 0f 1e fa          	endbr64 
    bafe:	55                   	push   %rbp
    baff:	48 89 e5             	mov    %rsp,%rbp
    bb02:	48 83 ec 10          	sub    $0x10,%rsp
    bb06:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bb0a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bb0e:	48 89 c7             	mov    %rax,%rdi
    bb11:	e8 c6 ff ff ff       	call   badc <_ZN3MPI5GroupD1Ev>
    bb16:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bb1a:	be 10 00 00 00       	mov    $0x10,%esi
    bb1f:	48 89 c7             	mov    %rax,%rdi
    bb22:	e8 c9 e8 ff ff       	call   a3f0 <_ZdlPvm@plt>
    bb27:	c9                   	leave  
    bb28:	c3                   	ret    
    bb29:	90                   	nop

000000000000bb2a <_ZNK3MPI5GroupcvP12ompi_group_tEv>:
    bb2a:	f3 0f 1e fa          	endbr64 
    bb2e:	55                   	push   %rbp
    bb2f:	48 89 e5             	mov    %rsp,%rbp
    bb32:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bb36:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bb3a:	48 8b 40 08          	mov    0x8(%rax),%rax
    bb3e:	5d                   	pop    %rbp
    bb3f:	c3                   	ret    

000000000000bb40 <_ZN3MPI9Comm_NullC1EP19ompi_communicator_t>:
    bb40:	f3 0f 1e fa          	endbr64 
    bb44:	55                   	push   %rbp
    bb45:	48 89 e5             	mov    %rsp,%rbp
    bb48:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bb4c:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    bb50:	48 8d 15 79 06 01 00 	lea    0x10679(%rip),%rdx        # 1c1d0 <_ZTVN3MPI9Comm_NullE+0x10>
    bb57:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bb5b:	48 89 10             	mov    %rdx,(%rax)
    bb5e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bb62:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bb66:	48 89 50 08          	mov    %rdx,0x8(%rax)
    bb6a:	90                   	nop
    bb6b:	5d                   	pop    %rbp
    bb6c:	c3                   	ret    
    bb6d:	90                   	nop

000000000000bb6e <_ZN3MPI9Comm_NullD1Ev>:
    bb6e:	f3 0f 1e fa          	endbr64 
    bb72:	55                   	push   %rbp
    bb73:	48 89 e5             	mov    %rsp,%rbp
    bb76:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bb7a:	48 8d 15 4f 06 01 00 	lea    0x1064f(%rip),%rdx        # 1c1d0 <_ZTVN3MPI9Comm_NullE+0x10>
    bb81:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bb85:	48 89 10             	mov    %rdx,(%rax)
    bb88:	90                   	nop
    bb89:	5d                   	pop    %rbp
    bb8a:	c3                   	ret    
    bb8b:	90                   	nop

000000000000bb8c <_ZN3MPI9Comm_NullD0Ev>:
    bb8c:	f3 0f 1e fa          	endbr64 
    bb90:	55                   	push   %rbp
    bb91:	48 89 e5             	mov    %rsp,%rbp
    bb94:	48 83 ec 10          	sub    $0x10,%rsp
    bb98:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bb9c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bba0:	48 89 c7             	mov    %rax,%rdi
    bba3:	e8 c6 ff ff ff       	call   bb6e <_ZN3MPI9Comm_NullD1Ev>
    bba8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bbac:	be 10 00 00 00       	mov    $0x10,%esi
    bbb1:	48 89 c7             	mov    %rax,%rdi
    bbb4:	e8 37 e8 ff ff       	call   a3f0 <_ZdlPvm@plt>
    bbb9:	c9                   	leave  
    bbba:	c3                   	ret    
    bbbb:	90                   	nop

000000000000bbbc <_ZNK3MPI9Comm_NullcvP19ompi_communicator_tEv>:
    bbbc:	f3 0f 1e fa          	endbr64 
    bbc0:	55                   	push   %rbp
    bbc1:	48 89 e5             	mov    %rsp,%rbp
    bbc4:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bbc8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bbcc:	48 8b 40 08          	mov    0x8(%rax),%rax
    bbd0:	5d                   	pop    %rbp
    bbd1:	c3                   	ret    

000000000000bbd2 <_ZN3MPI4CommC1EP19ompi_communicator_t>:
    bbd2:	f3 0f 1e fa          	endbr64 
    bbd6:	55                   	push   %rbp
    bbd7:	48 89 e5             	mov    %rsp,%rbp
    bbda:	48 83 ec 10          	sub    $0x10,%rsp
    bbde:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bbe2:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    bbe6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bbea:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bbee:	48 89 d6             	mov    %rdx,%rsi
    bbf1:	48 89 c7             	mov    %rax,%rdi
    bbf4:	e8 47 ff ff ff       	call   bb40 <_ZN3MPI9Comm_NullC1EP19ompi_communicator_t>
    bbf9:	48 8d 15 b0 06 01 00 	lea    0x106b0(%rip),%rdx        # 1c2b0 <_ZTVN3MPI4CommE+0x10>
    bc00:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bc04:	48 89 10             	mov    %rdx,(%rax)
    bc07:	90                   	nop
    bc08:	c9                   	leave  
    bc09:	c3                   	ret    

000000000000bc0a <_ZN3MPI3WinD1Ev>:
    bc0a:	f3 0f 1e fa          	endbr64 
    bc0e:	55                   	push   %rbp
    bc0f:	48 89 e5             	mov    %rsp,%rbp
    bc12:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bc16:	48 8d 15 d3 05 01 00 	lea    0x105d3(%rip),%rdx        # 1c1f0 <_ZTVN3MPI3WinE+0x10>
    bc1d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bc21:	48 89 10             	mov    %rdx,(%rax)
    bc24:	90                   	nop
    bc25:	5d                   	pop    %rbp
    bc26:	c3                   	ret    
    bc27:	90                   	nop

000000000000bc28 <_ZN3MPI3WinD0Ev>:
    bc28:	f3 0f 1e fa          	endbr64 
    bc2c:	55                   	push   %rbp
    bc2d:	48 89 e5             	mov    %rsp,%rbp
    bc30:	48 83 ec 10          	sub    $0x10,%rsp
    bc34:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bc38:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bc3c:	48 89 c7             	mov    %rax,%rdi
    bc3f:	e8 c6 ff ff ff       	call   bc0a <_ZN3MPI3WinD1Ev>
    bc44:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bc48:	be 10 00 00 00       	mov    $0x10,%esi
    bc4d:	48 89 c7             	mov    %rax,%rdi
    bc50:	e8 9b e7 ff ff       	call   a3f0 <_ZdlPvm@plt>
    bc55:	c9                   	leave  
    bc56:	c3                   	ret    
    bc57:	90                   	nop

000000000000bc58 <_ZN3MPI10ErrhandlerD1Ev>:
    bc58:	f3 0f 1e fa          	endbr64 
    bc5c:	55                   	push   %rbp
    bc5d:	48 89 e5             	mov    %rsp,%rbp
    bc60:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bc64:	48 8d 15 75 fa 00 00 	lea    0xfa75(%rip),%rdx        # 1b6e0 <_ZTVN3MPI10ErrhandlerE+0x10>
    bc6b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bc6f:	48 89 10             	mov    %rdx,(%rax)
    bc72:	90                   	nop
    bc73:	5d                   	pop    %rbp
    bc74:	c3                   	ret    
    bc75:	90                   	nop

000000000000bc76 <_ZN3MPI10ErrhandlerD0Ev>:
    bc76:	f3 0f 1e fa          	endbr64 
    bc7a:	55                   	push   %rbp
    bc7b:	48 89 e5             	mov    %rsp,%rbp
    bc7e:	48 83 ec 10          	sub    $0x10,%rsp
    bc82:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bc86:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bc8a:	48 89 c7             	mov    %rax,%rdi
    bc8d:	e8 c6 ff ff ff       	call   bc58 <_ZN3MPI10ErrhandlerD1Ev>
    bc92:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bc96:	be 10 00 00 00       	mov    $0x10,%esi
    bc9b:	48 89 c7             	mov    %rax,%rdi
    bc9e:	e8 4d e7 ff ff       	call   a3f0 <_ZdlPvm@plt>
    bca3:	c9                   	leave  
    bca4:	c3                   	ret    
    bca5:	90                   	nop

000000000000bca6 <_ZN3MPI10ErrhandlerC1EP17ompi_errhandler_t>:
    bca6:	f3 0f 1e fa          	endbr64 
    bcaa:	55                   	push   %rbp
    bcab:	48 89 e5             	mov    %rsp,%rbp
    bcae:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bcb2:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    bcb6:	48 8d 15 23 fa 00 00 	lea    0xfa23(%rip),%rdx        # 1b6e0 <_ZTVN3MPI10ErrhandlerE+0x10>
    bcbd:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bcc1:	48 89 10             	mov    %rdx,(%rax)
    bcc4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bcc8:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bccc:	48 89 50 08          	mov    %rdx,0x8(%rax)
    bcd0:	90                   	nop
    bcd1:	5d                   	pop    %rbp
    bcd2:	c3                   	ret    
    bcd3:	90                   	nop

000000000000bcd4 <_ZNK3MPI10ErrhandlercvP17ompi_errhandler_tEv>:
    bcd4:	f3 0f 1e fa          	endbr64 
    bcd8:	55                   	push   %rbp
    bcd9:	48 89 e5             	mov    %rsp,%rbp
    bcdc:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bce0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bce4:	48 8b 40 08          	mov    0x8(%rax),%rax
    bce8:	5d                   	pop    %rbp
    bce9:	c3                   	ret    

000000000000bcea <_ZN3MPI4CommD1Ev>:
    bcea:	f3 0f 1e fa          	endbr64 
    bcee:	55                   	push   %rbp
    bcef:	48 89 e5             	mov    %rsp,%rbp
    bcf2:	48 83 ec 10          	sub    $0x10,%rsp
    bcf6:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bcfa:	48 8d 15 af 05 01 00 	lea    0x105af(%rip),%rdx        # 1c2b0 <_ZTVN3MPI4CommE+0x10>
    bd01:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd05:	48 89 10             	mov    %rdx,(%rax)
    bd08:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd0c:	48 89 c7             	mov    %rax,%rdi
    bd0f:	e8 5a fe ff ff       	call   bb6e <_ZN3MPI9Comm_NullD1Ev>
    bd14:	90                   	nop
    bd15:	c9                   	leave  
    bd16:	c3                   	ret    
    bd17:	90                   	nop

000000000000bd18 <_ZN3MPI4CommD0Ev>:
    bd18:	f3 0f 1e fa          	endbr64 
    bd1c:	55                   	push   %rbp
    bd1d:	48 89 e5             	mov    %rsp,%rbp
    bd20:	48 83 ec 10          	sub    $0x10,%rsp
    bd24:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bd28:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd2c:	48 89 c7             	mov    %rax,%rdi
    bd2f:	e8 b6 ff ff ff       	call   bcea <_ZN3MPI4CommD1Ev>
    bd34:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd38:	be 10 00 00 00       	mov    $0x10,%esi
    bd3d:	48 89 c7             	mov    %rax,%rdi
    bd40:	e8 ab e6 ff ff       	call   a3f0 <_ZdlPvm@plt>
    bd45:	c9                   	leave  
    bd46:	c3                   	ret    
    bd47:	90                   	nop

000000000000bd48 <_ZN3MPI9IntracommC1Ev>:
    bd48:	f3 0f 1e fa          	endbr64 
    bd4c:	55                   	push   %rbp
    bd4d:	48 89 e5             	mov    %rsp,%rbp
    bd50:	48 83 ec 10          	sub    $0x10,%rsp
    bd54:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bd58:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd5c:	48 89 c7             	mov    %rax,%rdi
    bd5f:	e8 2c e3 ff ff       	call   a090 <_ZN3MPI4CommC2Ev@plt>
    bd64:	48 8d 15 fd 00 01 00 	lea    0x100fd(%rip),%rdx        # 1be68 <_ZTVN3MPI9IntracommE+0x10>
    bd6b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd6f:	48 89 10             	mov    %rdx,(%rax)
    bd72:	90                   	nop
    bd73:	c9                   	leave  
    bd74:	c3                   	ret    
    bd75:	90                   	nop

000000000000bd76 <_ZN3MPI9IntracommD1Ev>:
    bd76:	f3 0f 1e fa          	endbr64 
    bd7a:	55                   	push   %rbp
    bd7b:	48 89 e5             	mov    %rsp,%rbp
    bd7e:	48 83 ec 10          	sub    $0x10,%rsp
    bd82:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bd86:	48 8d 15 db 00 01 00 	lea    0x100db(%rip),%rdx        # 1be68 <_ZTVN3MPI9IntracommE+0x10>
    bd8d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd91:	48 89 10             	mov    %rdx,(%rax)
    bd94:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bd98:	48 89 c7             	mov    %rax,%rdi
    bd9b:	e8 4a ff ff ff       	call   bcea <_ZN3MPI4CommD1Ev>
    bda0:	90                   	nop
    bda1:	c9                   	leave  
    bda2:	c3                   	ret    
    bda3:	90                   	nop

000000000000bda4 <_ZN3MPI9IntracommD0Ev>:
    bda4:	f3 0f 1e fa          	endbr64 
    bda8:	55                   	push   %rbp
    bda9:	48 89 e5             	mov    %rsp,%rbp
    bdac:	48 83 ec 10          	sub    $0x10,%rsp
    bdb0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bdb4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bdb8:	48 89 c7             	mov    %rax,%rdi
    bdbb:	e8 b6 ff ff ff       	call   bd76 <_ZN3MPI9IntracommD1Ev>
    bdc0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bdc4:	be 10 00 00 00       	mov    $0x10,%esi
    bdc9:	48 89 c7             	mov    %rax,%rdi
    bdcc:	e8 1f e6 ff ff       	call   a3f0 <_ZdlPvm@plt>
    bdd1:	c9                   	leave  
    bdd2:	c3                   	ret    
    bdd3:	90                   	nop

000000000000bdd4 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>:
    bdd4:	f3 0f 1e fa          	endbr64 
    bdd8:	55                   	push   %rbp
    bdd9:	48 89 e5             	mov    %rsp,%rbp
    bddc:	48 83 ec 10          	sub    $0x10,%rsp
    bde0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    bde4:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    bde8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bdec:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bdf0:	48 89 d6             	mov    %rdx,%rsi
    bdf3:	48 89 c7             	mov    %rax,%rdi
    bdf6:	e8 d7 fd ff ff       	call   bbd2 <_ZN3MPI4CommC1EP19ompi_communicator_t>
    bdfb:	48 8d 15 9e f9 00 00 	lea    0xf99e(%rip),%rdx        # 1b7a0 <_ZTVN3MPI9IntercommE+0x10>
    be02:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    be06:	48 89 10             	mov    %rdx,(%rax)
    be09:	90                   	nop
    be0a:	c9                   	leave  
    be0b:	c3                   	ret    

000000000000be0c <_ZN3MPI4InfoC1EP11ompi_info_t>:
    be0c:	f3 0f 1e fa          	endbr64 
    be10:	55                   	push   %rbp
    be11:	48 89 e5             	mov    %rsp,%rbp
    be14:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    be18:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    be1c:	48 8d 15 e5 f7 00 00 	lea    0xf7e5(%rip),%rdx        # 1b608 <_ZTVN3MPI4InfoE+0x10>
    be23:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    be27:	48 89 10             	mov    %rdx,(%rax)
    be2a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    be2e:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    be32:	48 89 50 08          	mov    %rdx,0x8(%rax)
    be36:	90                   	nop
    be37:	5d                   	pop    %rbp
    be38:	c3                   	ret    
    be39:	90                   	nop

000000000000be3a <_ZN3MPI4InfoD1Ev>:
    be3a:	f3 0f 1e fa          	endbr64 
    be3e:	55                   	push   %rbp
    be3f:	48 89 e5             	mov    %rsp,%rbp
    be42:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    be46:	48 8d 15 bb f7 00 00 	lea    0xf7bb(%rip),%rdx        # 1b608 <_ZTVN3MPI4InfoE+0x10>
    be4d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    be51:	48 89 10             	mov    %rdx,(%rax)
    be54:	90                   	nop
    be55:	5d                   	pop    %rbp
    be56:	c3                   	ret    
    be57:	90                   	nop

000000000000be58 <_ZN3MPI4InfoD0Ev>:
    be58:	f3 0f 1e fa          	endbr64 
    be5c:	55                   	push   %rbp
    be5d:	48 89 e5             	mov    %rsp,%rbp
    be60:	48 83 ec 10          	sub    $0x10,%rsp
    be64:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    be68:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    be6c:	48 89 c7             	mov    %rax,%rdi
    be6f:	e8 c6 ff ff ff       	call   be3a <_ZN3MPI4InfoD1Ev>
    be74:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    be78:	be 10 00 00 00       	mov    $0x10,%esi
    be7d:	48 89 c7             	mov    %rax,%rdi
    be80:	e8 6b e5 ff ff       	call   a3f0 <_ZdlPvm@plt>
    be85:	c9                   	leave  
    be86:	c3                   	ret    
    be87:	90                   	nop

000000000000be88 <_ZNK3MPI4InfocvP11ompi_info_tEv>:
    be88:	f3 0f 1e fa          	endbr64 
    be8c:	55                   	push   %rbp
    be8d:	48 89 e5             	mov    %rsp,%rbp
    be90:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    be94:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    be98:	48 8b 40 08          	mov    0x8(%rax),%rax
    be9c:	5d                   	pop    %rbp
    be9d:	c3                   	ret    

000000000000be9e <_ZNK3MPI8Datatype17Create_contiguousEi>:
    be9e:	f3 0f 1e fa          	endbr64 
    bea2:	55                   	push   %rbp
    bea3:	48 89 e5             	mov    %rsp,%rbp
    bea6:	48 83 ec 30          	sub    $0x30,%rsp
    beaa:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    beae:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    beb2:	89 55 dc             	mov    %edx,-0x24(%rbp)
    beb5:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    bebc:	00 00 
    bebe:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    bec2:	31 c0                	xor    %eax,%eax
    bec4:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    bec8:	48 8b 48 08          	mov    0x8(%rax),%rcx
    becc:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    bed0:	8b 45 dc             	mov    -0x24(%rbp),%eax
    bed3:	48 89 ce             	mov    %rcx,%rsi
    bed6:	89 c7                	mov    %eax,%edi
    bed8:	e8 33 e2 ff ff       	call   a110 <MPI_Type_contiguous@plt>
    bedd:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bee1:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    bee5:	48 89 d6             	mov    %rdx,%rsi
    bee8:	48 89 c7             	mov    %rax,%rdi
    beeb:	e8 5c f9 ff ff       	call   b84c <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
    bef0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bef4:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    befb:	00 00 
    befd:	74 05                	je     bf04 <_ZNK3MPI8Datatype17Create_contiguousEi+0x66>
    beff:	e8 cc e5 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    bf04:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    bf08:	c9                   	leave  
    bf09:	c3                   	ret    

000000000000bf0a <_ZNK3MPI8Datatype13Create_vectorEiii>:
    bf0a:	f3 0f 1e fa          	endbr64 
    bf0e:	55                   	push   %rbp
    bf0f:	48 89 e5             	mov    %rsp,%rbp
    bf12:	48 83 ec 30          	sub    $0x30,%rsp
    bf16:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    bf1a:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    bf1e:	89 55 dc             	mov    %edx,-0x24(%rbp)
    bf21:	89 4d d8             	mov    %ecx,-0x28(%rbp)
    bf24:	44 89 45 d4          	mov    %r8d,-0x2c(%rbp)
    bf28:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    bf2f:	00 00 
    bf31:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    bf35:	31 c0                	xor    %eax,%eax
    bf37:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    bf3b:	48 8b 48 08          	mov    0x8(%rax),%rcx
    bf3f:	48 8d 7d f0          	lea    -0x10(%rbp),%rdi
    bf43:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    bf46:	8b 75 d8             	mov    -0x28(%rbp),%esi
    bf49:	8b 45 dc             	mov    -0x24(%rbp),%eax
    bf4c:	49 89 f8             	mov    %rdi,%r8
    bf4f:	89 c7                	mov    %eax,%edi
    bf51:	e8 8a e9 ff ff       	call   a8e0 <MPI_Type_vector@plt>
    bf56:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bf5a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    bf5e:	48 89 d6             	mov    %rdx,%rsi
    bf61:	48 89 c7             	mov    %rax,%rdi
    bf64:	e8 e3 f8 ff ff       	call   b84c <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
    bf69:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bf6d:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    bf74:	00 00 
    bf76:	74 05                	je     bf7d <_ZNK3MPI8Datatype13Create_vectorEiii+0x73>
    bf78:	e8 53 e5 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    bf7d:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    bf81:	c9                   	leave  
    bf82:	c3                   	ret    
    bf83:	90                   	nop

000000000000bf84 <_ZNK3MPI8Datatype14Create_indexedEiPKiS2_>:
    bf84:	f3 0f 1e fa          	endbr64 
    bf88:	55                   	push   %rbp
    bf89:	48 89 e5             	mov    %rsp,%rbp
    bf8c:	48 83 ec 40          	sub    $0x40,%rsp
    bf90:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    bf94:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    bf98:	89 55 dc             	mov    %edx,-0x24(%rbp)
    bf9b:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    bf9f:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    bfa3:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    bfaa:	00 00 
    bfac:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    bfb0:	31 c0                	xor    %eax,%eax
    bfb2:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    bfb6:	48 8b 48 08          	mov    0x8(%rax),%rcx
    bfba:	48 8d 7d f0          	lea    -0x10(%rbp),%rdi
    bfbe:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    bfc2:	48 8b 75 d0          	mov    -0x30(%rbp),%rsi
    bfc6:	8b 45 dc             	mov    -0x24(%rbp),%eax
    bfc9:	49 89 f8             	mov    %rdi,%r8
    bfcc:	89 c7                	mov    %eax,%edi
    bfce:	e8 3d e9 ff ff       	call   a910 <MPI_Type_indexed@plt>
    bfd3:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    bfd7:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    bfdb:	48 89 d6             	mov    %rdx,%rsi
    bfde:	48 89 c7             	mov    %rax,%rdi
    bfe1:	e8 66 f8 ff ff       	call   b84c <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
    bfe6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    bfea:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    bff1:	00 00 
    bff3:	74 05                	je     bffa <_ZNK3MPI8Datatype14Create_indexedEiPKiS2_+0x76>
    bff5:	e8 d6 e4 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    bffa:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    bffe:	c9                   	leave  
    bfff:	c3                   	ret    

000000000000c000 <_ZNK3MPI8Datatype15Create_hindexedEiPKiPKl>:
    c000:	f3 0f 1e fa          	endbr64 
    c004:	55                   	push   %rbp
    c005:	48 89 e5             	mov    %rsp,%rbp
    c008:	48 83 ec 40          	sub    $0x40,%rsp
    c00c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c010:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c014:	89 55 dc             	mov    %edx,-0x24(%rbp)
    c017:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    c01b:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    c01f:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c026:	00 00 
    c028:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c02c:	31 c0                	xor    %eax,%eax
    c02e:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c032:	48 8b 48 08          	mov    0x8(%rax),%rcx
    c036:	48 8d 7d f0          	lea    -0x10(%rbp),%rdi
    c03a:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    c03e:	48 8b 75 d0          	mov    -0x30(%rbp),%rsi
    c042:	8b 45 dc             	mov    -0x24(%rbp),%eax
    c045:	49 89 f8             	mov    %rdi,%r8
    c048:	89 c7                	mov    %eax,%edi
    c04a:	e8 b1 e4 ff ff       	call   a500 <MPI_Type_create_hindexed@plt>
    c04f:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c053:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c057:	48 89 d6             	mov    %rdx,%rsi
    c05a:	48 89 c7             	mov    %rax,%rdi
    c05d:	e8 ea f7 ff ff       	call   b84c <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
    c062:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c066:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c06d:	00 00 
    c06f:	74 05                	je     c076 <_ZNK3MPI8Datatype15Create_hindexedEiPKiPKl+0x76>
    c071:	e8 5a e4 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    c076:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c07a:	c9                   	leave  
    c07b:	c3                   	ret    

000000000000c07c <_ZNK3MPI8Datatype14Create_hvectorEiil>:
    c07c:	f3 0f 1e fa          	endbr64 
    c080:	55                   	push   %rbp
    c081:	48 89 e5             	mov    %rsp,%rbp
    c084:	48 83 ec 30          	sub    $0x30,%rsp
    c088:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c08c:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c090:	89 55 dc             	mov    %edx,-0x24(%rbp)
    c093:	89 4d d8             	mov    %ecx,-0x28(%rbp)
    c096:	4c 89 45 d0          	mov    %r8,-0x30(%rbp)
    c09a:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c0a1:	00 00 
    c0a3:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c0a7:	31 c0                	xor    %eax,%eax
    c0a9:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c0ad:	48 8b 48 08          	mov    0x8(%rax),%rcx
    c0b1:	48 8d 7d f0          	lea    -0x10(%rbp),%rdi
    c0b5:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    c0b9:	8b 75 d8             	mov    -0x28(%rbp),%esi
    c0bc:	8b 45 dc             	mov    -0x24(%rbp),%eax
    c0bf:	49 89 f8             	mov    %rdi,%r8
    c0c2:	89 c7                	mov    %eax,%edi
    c0c4:	e8 77 de ff ff       	call   9f40 <MPI_Type_create_hvector@plt>
    c0c9:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c0cd:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c0d1:	48 89 d6             	mov    %rdx,%rsi
    c0d4:	48 89 c7             	mov    %rax,%rdi
    c0d7:	e8 70 f7 ff ff       	call   b84c <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
    c0dc:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c0e0:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c0e7:	00 00 
    c0e9:	74 05                	je     c0f0 <_ZNK3MPI8Datatype14Create_hvectorEiil+0x74>
    c0eb:	e8 e0 e3 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    c0f0:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c0f4:	c9                   	leave  
    c0f5:	c3                   	ret    

000000000000c0f6 <_ZNK3MPI8Datatype20Create_indexed_blockEiiPKi>:
    c0f6:	f3 0f 1e fa          	endbr64 
    c0fa:	55                   	push   %rbp
    c0fb:	48 89 e5             	mov    %rsp,%rbp
    c0fe:	48 83 ec 30          	sub    $0x30,%rsp
    c102:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c106:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c10a:	89 55 dc             	mov    %edx,-0x24(%rbp)
    c10d:	89 4d d8             	mov    %ecx,-0x28(%rbp)
    c110:	4c 89 45 d0          	mov    %r8,-0x30(%rbp)
    c114:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c11b:	00 00 
    c11d:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c121:	31 c0                	xor    %eax,%eax
    c123:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c127:	48 8b 48 08          	mov    0x8(%rax),%rcx
    c12b:	48 8d 7d f0          	lea    -0x10(%rbp),%rdi
    c12f:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    c133:	8b 75 d8             	mov    -0x28(%rbp),%esi
    c136:	8b 45 dc             	mov    -0x24(%rbp),%eax
    c139:	49 89 f8             	mov    %rdi,%r8
    c13c:	89 c7                	mov    %eax,%edi
    c13e:	e8 6d e5 ff ff       	call   a6b0 <MPI_Type_create_indexed_block@plt>
    c143:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c147:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c14b:	48 89 d6             	mov    %rdx,%rsi
    c14e:	48 89 c7             	mov    %rax,%rdi
    c151:	e8 f6 f6 ff ff       	call   b84c <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
    c156:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c15a:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c161:	00 00 
    c163:	74 05                	je     c16a <_ZNK3MPI8Datatype20Create_indexed_blockEiiPKi+0x74>
    c165:	e8 66 e3 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    c16a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c16e:	c9                   	leave  
    c16f:	c3                   	ret    

000000000000c170 <_ZNK3MPI8Datatype14Create_resizedEll>:
    c170:	f3 0f 1e fa          	endbr64 
    c174:	55                   	push   %rbp
    c175:	48 89 e5             	mov    %rsp,%rbp
    c178:	48 83 ec 30          	sub    $0x30,%rsp
    c17c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c180:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c184:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    c188:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    c18c:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c193:	00 00 
    c195:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c199:	31 c0                	xor    %eax,%eax
    c19b:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c19f:	48 8b 40 08          	mov    0x8(%rax),%rax
    c1a3:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    c1a7:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    c1ab:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    c1af:	48 89 c7             	mov    %rax,%rdi
    c1b2:	e8 b9 e7 ff ff       	call   a970 <MPI_Type_create_resized@plt>
    c1b7:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c1bb:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c1bf:	48 89 d6             	mov    %rdx,%rsi
    c1c2:	48 89 c7             	mov    %rax,%rdi
    c1c5:	e8 82 f6 ff ff       	call   b84c <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
    c1ca:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c1ce:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c1d5:	00 00 
    c1d7:	74 05                	je     c1de <_ZNK3MPI8Datatype14Create_resizedEll+0x6e>
    c1d9:	e8 f2 e2 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    c1de:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c1e2:	c9                   	leave  
    c1e3:	c3                   	ret    

000000000000c1e4 <_ZNK3MPI8Datatype8Get_sizeEv>:
    c1e4:	f3 0f 1e fa          	endbr64 
    c1e8:	55                   	push   %rbp
    c1e9:	48 89 e5             	mov    %rsp,%rbp
    c1ec:	48 83 ec 20          	sub    $0x20,%rsp
    c1f0:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c1f4:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c1fb:	00 00 
    c1fd:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c201:	31 c0                	xor    %eax,%eax
    c203:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c207:	48 8b 40 08          	mov    0x8(%rax),%rax
    c20b:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    c20f:	48 89 d6             	mov    %rdx,%rsi
    c212:	48 89 c7             	mov    %rax,%rdi
    c215:	e8 b6 e3 ff ff       	call   a5d0 <MPI_Type_size@plt>
    c21a:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c21d:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    c221:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    c228:	00 00 
    c22a:	74 05                	je     c231 <_ZNK3MPI8Datatype8Get_sizeEv+0x4d>
    c22c:	e8 9f e2 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    c231:	c9                   	leave  
    c232:	c3                   	ret    
    c233:	90                   	nop

000000000000c234 <_ZNK3MPI8Datatype10Get_extentERlS1_>:
    c234:	f3 0f 1e fa          	endbr64 
    c238:	55                   	push   %rbp
    c239:	48 89 e5             	mov    %rsp,%rbp
    c23c:	48 83 ec 20          	sub    $0x20,%rsp
    c240:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c244:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c248:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    c24c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c250:	48 8b 40 08          	mov    0x8(%rax),%rax
    c254:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    c258:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    c25c:	48 89 ce             	mov    %rcx,%rsi
    c25f:	48 89 c7             	mov    %rax,%rdi
    c262:	e8 69 dd ff ff       	call   9fd0 <MPI_Type_get_extent@plt>
    c267:	90                   	nop
    c268:	c9                   	leave  
    c269:	c3                   	ret    

000000000000c26a <_ZNK3MPI8Datatype15Get_true_extentERlS1_>:
    c26a:	f3 0f 1e fa          	endbr64 
    c26e:	55                   	push   %rbp
    c26f:	48 89 e5             	mov    %rsp,%rbp
    c272:	48 83 ec 20          	sub    $0x20,%rsp
    c276:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c27a:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c27e:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    c282:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c286:	48 8b 40 08          	mov    0x8(%rax),%rax
    c28a:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    c28e:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    c292:	48 89 ce             	mov    %rcx,%rsi
    c295:	48 89 c7             	mov    %rax,%rdi
    c298:	e8 43 de ff ff       	call   a0e0 <MPI_Type_get_true_extent@plt>
    c29d:	90                   	nop
    c29e:	c9                   	leave  
    c29f:	c3                   	ret    

000000000000c2a0 <_ZN3MPI8Datatype6CommitEv>:
    c2a0:	f3 0f 1e fa          	endbr64 
    c2a4:	55                   	push   %rbp
    c2a5:	48 89 e5             	mov    %rsp,%rbp
    c2a8:	48 83 ec 10          	sub    $0x10,%rsp
    c2ac:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c2b0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c2b4:	48 83 c0 08          	add    $0x8,%rax
    c2b8:	48 89 c7             	mov    %rax,%rdi
    c2bb:	e8 e0 e2 ff ff       	call   a5a0 <MPI_Type_commit@plt>
    c2c0:	90                   	nop
    c2c1:	c9                   	leave  
    c2c2:	c3                   	ret    
    c2c3:	90                   	nop

000000000000c2c4 <_ZNK3MPI8Datatype4PackEPKviPviRiRKNS_4CommE>:
    c2c4:	f3 0f 1e fa          	endbr64 
    c2c8:	55                   	push   %rbp
    c2c9:	48 89 e5             	mov    %rsp,%rbp
    c2cc:	48 83 ec 30          	sub    $0x30,%rsp
    c2d0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c2d4:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c2d8:	89 55 ec             	mov    %edx,-0x14(%rbp)
    c2db:	48 89 4d e0          	mov    %rcx,-0x20(%rbp)
    c2df:	44 89 45 e8          	mov    %r8d,-0x18(%rbp)
    c2e3:	4c 89 4d d8          	mov    %r9,-0x28(%rbp)
    c2e7:	48 8b 45 10          	mov    0x10(%rbp),%rax
    c2eb:	48 89 c7             	mov    %rax,%rdi
    c2ee:	e8 c9 f8 ff ff       	call   bbbc <_ZNK3MPI9Comm_NullcvP19ompi_communicator_tEv>
    c2f3:	48 89 c7             	mov    %rax,%rdi
    c2f6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c2fa:	48 8b 50 08          	mov    0x8(%rax),%rdx
    c2fe:	4c 8b 4d d8          	mov    -0x28(%rbp),%r9
    c302:	44 8b 45 e8          	mov    -0x18(%rbp),%r8d
    c306:	48 8b 4d e0          	mov    -0x20(%rbp),%rcx
    c30a:	8b 75 ec             	mov    -0x14(%rbp),%esi
    c30d:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    c311:	48 83 ec 08          	sub    $0x8,%rsp
    c315:	57                   	push   %rdi
    c316:	48 89 c7             	mov    %rax,%rdi
    c319:	e8 82 e0 ff ff       	call   a3a0 <MPI_Pack@plt>
    c31e:	48 83 c4 10          	add    $0x10,%rsp
    c322:	90                   	nop
    c323:	c9                   	leave  
    c324:	c3                   	ret    
    c325:	90                   	nop

000000000000c326 <_ZNK3MPI8Datatype6UnpackEPKviPviRiRKNS_4CommE>:
    c326:	f3 0f 1e fa          	endbr64 
    c32a:	55                   	push   %rbp
    c32b:	48 89 e5             	mov    %rsp,%rbp
    c32e:	48 83 ec 30          	sub    $0x30,%rsp
    c332:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c336:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c33a:	89 55 ec             	mov    %edx,-0x14(%rbp)
    c33d:	48 89 4d e0          	mov    %rcx,-0x20(%rbp)
    c341:	44 89 45 e8          	mov    %r8d,-0x18(%rbp)
    c345:	4c 89 4d d8          	mov    %r9,-0x28(%rbp)
    c349:	48 8b 45 10          	mov    0x10(%rbp),%rax
    c34d:	48 89 c7             	mov    %rax,%rdi
    c350:	e8 67 f8 ff ff       	call   bbbc <_ZNK3MPI9Comm_NullcvP19ompi_communicator_tEv>
    c355:	48 89 c7             	mov    %rax,%rdi
    c358:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c35c:	4c 8b 48 08          	mov    0x8(%rax),%r9
    c360:	44 8b 45 e8          	mov    -0x18(%rbp),%r8d
    c364:	48 8b 4d e0          	mov    -0x20(%rbp),%rcx
    c368:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
    c36c:	8b 75 ec             	mov    -0x14(%rbp),%esi
    c36f:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    c373:	48 83 ec 08          	sub    $0x8,%rsp
    c377:	57                   	push   %rdi
    c378:	48 89 c7             	mov    %rax,%rdi
    c37b:	e8 80 dc ff ff       	call   a000 <MPI_Unpack@plt>
    c380:	48 83 c4 10          	add    $0x10,%rsp
    c384:	90                   	nop
    c385:	c9                   	leave  
    c386:	c3                   	ret    
    c387:	90                   	nop

000000000000c388 <_ZNK3MPI8Datatype9Pack_sizeEiRKNS_4CommE>:
    c388:	f3 0f 1e fa          	endbr64 
    c38c:	55                   	push   %rbp
    c38d:	48 89 e5             	mov    %rsp,%rbp
    c390:	48 83 ec 30          	sub    $0x30,%rsp
    c394:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c398:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    c39b:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    c39f:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c3a6:	00 00 
    c3a8:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c3ac:	31 c0                	xor    %eax,%eax
    c3ae:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    c3b2:	48 89 c7             	mov    %rax,%rdi
    c3b5:	e8 02 f8 ff ff       	call   bbbc <_ZNK3MPI9Comm_NullcvP19ompi_communicator_tEv>
    c3ba:	48 89 c7             	mov    %rax,%rdi
    c3bd:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c3c1:	48 8b 70 08          	mov    0x8(%rax),%rsi
    c3c5:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    c3c9:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    c3cc:	48 89 d1             	mov    %rdx,%rcx
    c3cf:	48 89 fa             	mov    %rdi,%rdx
    c3d2:	89 c7                	mov    %eax,%edi
    c3d4:	e8 37 db ff ff       	call   9f10 <MPI_Pack_size@plt>
    c3d9:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c3dc:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    c3e0:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    c3e7:	00 00 
    c3e9:	74 05                	je     c3f0 <_ZNK3MPI8Datatype9Pack_sizeEiRKNS_4CommE+0x68>
    c3eb:	e8 e0 e0 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    c3f0:	c9                   	leave  
    c3f1:	c3                   	ret    

000000000000c3f2 <_ZNK3MPI8Datatype13Pack_externalEPKcPKviPvlRl>:
    c3f2:	f3 0f 1e fa          	endbr64 
    c3f6:	55                   	push   %rbp
    c3f7:	48 89 e5             	mov    %rsp,%rbp
    c3fa:	48 83 ec 30          	sub    $0x30,%rsp
    c3fe:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c402:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c406:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    c40a:	89 4d e4             	mov    %ecx,-0x1c(%rbp)
    c40d:	4c 89 45 d8          	mov    %r8,-0x28(%rbp)
    c411:	4c 89 4d d0          	mov    %r9,-0x30(%rbp)
    c415:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c419:	48 8b 48 08          	mov    0x8(%rax),%rcx
    c41d:	4c 8b 45 d0          	mov    -0x30(%rbp),%r8
    c421:	48 8b 7d d8          	mov    -0x28(%rbp),%rdi
    c425:	8b 55 e4             	mov    -0x1c(%rbp),%edx
    c428:	48 8b 75 e8          	mov    -0x18(%rbp),%rsi
    c42c:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    c430:	48 83 ec 08          	sub    $0x8,%rsp
    c434:	ff 75 10             	push   0x10(%rbp)
    c437:	4d 89 c1             	mov    %r8,%r9
    c43a:	49 89 f8             	mov    %rdi,%r8
    c43d:	48 89 c7             	mov    %rax,%rdi
    c440:	e8 9b dd ff ff       	call   a1e0 <MPI_Pack_external@plt>
    c445:	48 83 c4 10          	add    $0x10,%rsp
    c449:	90                   	nop
    c44a:	c9                   	leave  
    c44b:	c3                   	ret    

000000000000c44c <_ZNK3MPI8Datatype18Pack_external_sizeEPKci>:
    c44c:	f3 0f 1e fa          	endbr64 
    c450:	55                   	push   %rbp
    c451:	48 89 e5             	mov    %rsp,%rbp
    c454:	48 83 ec 30          	sub    $0x30,%rsp
    c458:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c45c:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c460:	89 55 dc             	mov    %edx,-0x24(%rbp)
    c463:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c46a:	00 00 
    c46c:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c470:	31 c0                	xor    %eax,%eax
    c472:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c476:	48 8b 50 08          	mov    0x8(%rax),%rdx
    c47a:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    c47e:	8b 75 dc             	mov    -0x24(%rbp),%esi
    c481:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c485:	48 89 c7             	mov    %rax,%rdi
    c488:	e8 03 e6 ff ff       	call   aa90 <MPI_Pack_external_size@plt>
    c48d:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    c491:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    c495:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    c49c:	00 00 
    c49e:	74 05                	je     c4a5 <_ZNK3MPI8Datatype18Pack_external_sizeEPKci+0x59>
    c4a0:	e8 2b e0 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    c4a5:	c9                   	leave  
    c4a6:	c3                   	ret    
    c4a7:	90                   	nop

000000000000c4a8 <_ZNK3MPI8Datatype15Unpack_externalEPKcPKvlRlPvi>:
    c4a8:	f3 0f 1e fa          	endbr64 
    c4ac:	55                   	push   %rbp
    c4ad:	48 89 e5             	mov    %rsp,%rbp
    c4b0:	48 83 ec 30          	sub    $0x30,%rsp
    c4b4:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c4b8:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c4bc:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    c4c0:	48 89 4d e0          	mov    %rcx,-0x20(%rbp)
    c4c4:	4c 89 45 d8          	mov    %r8,-0x28(%rbp)
    c4c8:	4c 89 4d d0          	mov    %r9,-0x30(%rbp)
    c4cc:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c4d0:	48 8b 78 08          	mov    0x8(%rax),%rdi
    c4d4:	4c 8b 45 d0          	mov    -0x30(%rbp),%r8
    c4d8:	48 8b 4d d8          	mov    -0x28(%rbp),%rcx
    c4dc:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    c4e0:	48 8b 75 e8          	mov    -0x18(%rbp),%rsi
    c4e4:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    c4e8:	48 83 ec 08          	sub    $0x8,%rsp
    c4ec:	57                   	push   %rdi
    c4ed:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
    c4f1:	48 89 c7             	mov    %rax,%rdi
    c4f4:	e8 57 da ff ff       	call   9f50 <MPI_Unpack_external@plt>
    c4f9:	48 83 c4 10          	add    $0x10,%rsp
    c4fd:	90                   	nop
    c4fe:	c9                   	leave  
    c4ff:	c3                   	ret    

000000000000c500 <_ZNK3MPI8Datatype15Create_subarrayEiPKiS2_S2_i>:
    c500:	f3 0f 1e fa          	endbr64 
    c504:	55                   	push   %rbp
    c505:	48 89 e5             	mov    %rsp,%rbp
    c508:	48 83 ec 40          	sub    $0x40,%rsp
    c50c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c510:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c514:	89 55 dc             	mov    %edx,-0x24(%rbp)
    c517:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    c51b:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    c51f:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    c523:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c52a:	00 00 
    c52c:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c530:	31 c0                	xor    %eax,%eax
    c532:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c536:	4c 8b 40 08          	mov    0x8(%rax),%r8
    c53a:	48 8b 4d c0          	mov    -0x40(%rbp),%rcx
    c53e:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    c542:	48 8b 75 d0          	mov    -0x30(%rbp),%rsi
    c546:	8b 45 dc             	mov    -0x24(%rbp),%eax
    c549:	48 83 ec 08          	sub    $0x8,%rsp
    c54d:	48 8d 7d f0          	lea    -0x10(%rbp),%rdi
    c551:	57                   	push   %rdi
    c552:	4d 89 c1             	mov    %r8,%r9
    c555:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    c559:	89 c7                	mov    %eax,%edi
    c55b:	e8 20 e3 ff ff       	call   a880 <MPI_Type_create_subarray@plt>
    c560:	48 83 c4 10          	add    $0x10,%rsp
    c564:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c568:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c56c:	48 89 d6             	mov    %rdx,%rsi
    c56f:	48 89 c7             	mov    %rax,%rdi
    c572:	e8 d5 f2 ff ff       	call   b84c <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
    c577:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c57b:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c582:	00 00 
    c584:	74 05                	je     c58b <_ZNK3MPI8Datatype15Create_subarrayEiPKiS2_S2_i+0x8b>
    c586:	e8 45 df ff ff       	call   a4d0 <__stack_chk_fail@plt>
    c58b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c58f:	c9                   	leave  
    c590:	c3                   	ret    
    c591:	90                   	nop

000000000000c592 <_ZNK3MPI8Datatype13Create_darrayEiiiPKiS2_S2_S2_i>:
    c592:	f3 0f 1e fa          	endbr64 
    c596:	55                   	push   %rbp
    c597:	48 89 e5             	mov    %rsp,%rbp
    c59a:	48 83 ec 50          	sub    $0x50,%rsp
    c59e:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c5a2:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c5a6:	89 55 dc             	mov    %edx,-0x24(%rbp)
    c5a9:	89 4d d8             	mov    %ecx,-0x28(%rbp)
    c5ac:	44 89 45 d4          	mov    %r8d,-0x2c(%rbp)
    c5b0:	4c 89 4d c8          	mov    %r9,-0x38(%rbp)
    c5b4:	48 8b 45 10          	mov    0x10(%rbp),%rax
    c5b8:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
    c5bc:	48 8b 45 18          	mov    0x18(%rbp),%rax
    c5c0:	48 89 45 b8          	mov    %rax,-0x48(%rbp)
    c5c4:	48 8b 45 20          	mov    0x20(%rbp),%rax
    c5c8:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
    c5cc:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c5d3:	00 00 
    c5d5:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c5d9:	31 c0                	xor    %eax,%eax
    c5db:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c5df:	48 8b 78 08          	mov    0x8(%rax),%rdi
    c5e3:	4c 8b 4d b8          	mov    -0x48(%rbp),%r9
    c5e7:	4c 8b 55 c0          	mov    -0x40(%rbp),%r10
    c5eb:	48 8b 4d c8          	mov    -0x38(%rbp),%rcx
    c5ef:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    c5f2:	8b 75 d8             	mov    -0x28(%rbp),%esi
    c5f5:	8b 45 dc             	mov    -0x24(%rbp),%eax
    c5f8:	4c 8d 45 f0          	lea    -0x10(%rbp),%r8
    c5fc:	41 50                	push   %r8
    c5fe:	57                   	push   %rdi
    c5ff:	8b 7d 28             	mov    0x28(%rbp),%edi
    c602:	57                   	push   %rdi
    c603:	ff 75 b0             	push   -0x50(%rbp)
    c606:	4d 89 d0             	mov    %r10,%r8
    c609:	89 c7                	mov    %eax,%edi
    c60b:	e8 90 dc ff ff       	call   a2a0 <MPI_Type_create_darray@plt>
    c610:	48 83 c4 20          	add    $0x20,%rsp
    c614:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c618:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c61c:	48 89 d6             	mov    %rdx,%rsi
    c61f:	48 89 c7             	mov    %rax,%rdi
    c622:	e8 25 f2 ff ff       	call   b84c <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
    c627:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c62b:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c632:	00 00 
    c634:	74 05                	je     c63b <_ZNK3MPI8Datatype13Create_darrayEiiiPKiS2_S2_S2_i+0xa9>
    c636:	e8 95 de ff ff       	call   a4d0 <__stack_chk_fail@plt>
    c63b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c63f:	c9                   	leave  
    c640:	c3                   	ret    
    c641:	90                   	nop

000000000000c642 <_ZNK3MPI8Datatype3DupEv>:
    c642:	f3 0f 1e fa          	endbr64 
    c646:	55                   	push   %rbp
    c647:	48 89 e5             	mov    %rsp,%rbp
    c64a:	48 83 ec 20          	sub    $0x20,%rsp
    c64e:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c652:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c656:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c65d:	00 00 
    c65f:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c663:	31 c0                	xor    %eax,%eax
    c665:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c669:	48 8b 40 08          	mov    0x8(%rax),%rax
    c66d:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    c671:	48 89 d6             	mov    %rdx,%rsi
    c674:	48 89 c7             	mov    %rax,%rdi
    c677:	e8 24 e1 ff ff       	call   a7a0 <MPI_Type_dup@plt>
    c67c:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c680:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c684:	48 89 d6             	mov    %rdx,%rsi
    c687:	48 89 c7             	mov    %rax,%rdi
    c68a:	e8 bd f1 ff ff       	call   b84c <_ZN3MPI8DatatypeC1EP15ompi_datatype_t>
    c68f:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c693:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    c69a:	00 00 
    c69c:	74 05                	je     c6a3 <_ZNK3MPI8Datatype3DupEv+0x61>
    c69e:	e8 2d de ff ff       	call   a4d0 <__stack_chk_fail@plt>
    c6a3:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c6a7:	c9                   	leave  
    c6a8:	c3                   	ret    
    c6a9:	90                   	nop

000000000000c6aa <_ZN3MPI8Datatype11Delete_attrEi>:
    c6aa:	f3 0f 1e fa          	endbr64 
    c6ae:	55                   	push   %rbp
    c6af:	48 89 e5             	mov    %rsp,%rbp
    c6b2:	48 83 ec 10          	sub    $0x10,%rsp
    c6b6:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c6ba:	89 75 f4             	mov    %esi,-0xc(%rbp)
    c6bd:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c6c1:	48 8b 40 08          	mov    0x8(%rax),%rax
    c6c5:	8b 55 f4             	mov    -0xc(%rbp),%edx
    c6c8:	89 d6                	mov    %edx,%esi
    c6ca:	48 89 c7             	mov    %rax,%rdi
    c6cd:	e8 0e e3 ff ff       	call   a9e0 <MPI_Type_delete_attr@plt>
    c6d2:	90                   	nop
    c6d3:	c9                   	leave  
    c6d4:	c3                   	ret    
    c6d5:	90                   	nop

000000000000c6d6 <_ZNK3MPI8Datatype8Get_attrEiPv>:
    c6d6:	f3 0f 1e fa          	endbr64 
    c6da:	55                   	push   %rbp
    c6db:	48 89 e5             	mov    %rsp,%rbp
    c6de:	48 83 ec 30          	sub    $0x30,%rsp
    c6e2:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c6e6:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    c6e9:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    c6ed:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c6f4:	00 00 
    c6f6:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c6fa:	31 c0                	xor    %eax,%eax
    c6fc:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c700:	48 8b 40 08          	mov    0x8(%rax),%rax
    c704:	48 8d 4d f4          	lea    -0xc(%rbp),%rcx
    c708:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
    c70c:	8b 75 e4             	mov    -0x1c(%rbp),%esi
    c70f:	48 89 c7             	mov    %rax,%rdi
    c712:	e8 19 dc ff ff       	call   a330 <MPI_Type_get_attr@plt>
    c717:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c71a:	85 c0                	test   %eax,%eax
    c71c:	0f 95 c0             	setne  %al
    c71f:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    c723:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    c72a:	00 00 
    c72c:	74 05                	je     c733 <_ZNK3MPI8Datatype8Get_attrEiPv+0x5d>
    c72e:	e8 9d dd ff ff       	call   a4d0 <__stack_chk_fail@plt>
    c733:	c9                   	leave  
    c734:	c3                   	ret    
    c735:	90                   	nop

000000000000c736 <_ZNK3MPI8Datatype12Get_contentsEiiiPiPlPS0_>:
    c736:	f3 0f 1e fa          	endbr64 
    c73a:	55                   	push   %rbp
    c73b:	48 89 e5             	mov    %rsp,%rbp
    c73e:	48 83 ec 40          	sub    $0x40,%rsp
    c742:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c746:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    c749:	89 55 e0             	mov    %edx,-0x20(%rbp)
    c74c:	89 4d dc             	mov    %ecx,-0x24(%rbp)
    c74f:	4c 89 45 d0          	mov    %r8,-0x30(%rbp)
    c753:	4c 89 4d c8          	mov    %r9,-0x38(%rbp)
    c757:	8b 45 dc             	mov    -0x24(%rbp),%eax
    c75a:	48 98                	cltq   
    c75c:	48 ba ff ff ff ff ff 	movabs $0xfffffffffffffff,%rdx
    c763:	ff ff 0f 
    c766:	48 39 d0             	cmp    %rdx,%rax
    c769:	77 4b                	ja     c7b6 <_ZNK3MPI8Datatype12Get_contentsEiiiPiPlPS0_+0x80>
    c76b:	48 c1 e0 03          	shl    $0x3,%rax
    c76f:	48 89 c7             	mov    %rax,%rdi
    c772:	e8 29 d6 ff ff       	call   9da0 <_Znam@plt>
    c777:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c77b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c77f:	48 8b 40 08          	mov    0x8(%rax),%rax
    c783:	4c 8b 45 c8          	mov    -0x38(%rbp),%r8
    c787:	48 8b 7d d0          	mov    -0x30(%rbp),%rdi
    c78b:	8b 4d dc             	mov    -0x24(%rbp),%ecx
    c78e:	8b 55 e0             	mov    -0x20(%rbp),%edx
    c791:	8b 75 e4             	mov    -0x1c(%rbp),%esi
    c794:	48 83 ec 08          	sub    $0x8,%rsp
    c798:	ff 75 f8             	push   -0x8(%rbp)
    c79b:	4d 89 c1             	mov    %r8,%r9
    c79e:	49 89 f8             	mov    %rdi,%r8
    c7a1:	48 89 c7             	mov    %rax,%rdi
    c7a4:	e8 17 df ff ff       	call   a6c0 <MPI_Type_get_contents@plt>
    c7a9:	48 83 c4 10          	add    $0x10,%rsp
    c7ad:	c7 45 f4 00 00 00 00 	movl   $0x0,-0xc(%rbp)
    c7b4:	eb 3b                	jmp    c7f1 <_ZNK3MPI8Datatype12Get_contentsEiiiPiPlPS0_+0xbb>
    c7b6:	e8 25 df ff ff       	call   a6e0 <__cxa_throw_bad_array_new_length@plt>
    c7bb:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c7be:	48 98                	cltq   
    c7c0:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
    c7c7:	00 
    c7c8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c7cc:	48 01 c2             	add    %rax,%rdx
    c7cf:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c7d2:	48 98                	cltq   
    c7d4:	48 c1 e0 04          	shl    $0x4,%rax
    c7d8:	48 89 c1             	mov    %rax,%rcx
    c7db:	48 8b 45 10          	mov    0x10(%rbp),%rax
    c7df:	48 01 c8             	add    %rcx,%rax
    c7e2:	48 89 d6             	mov    %rdx,%rsi
    c7e5:	48 89 c7             	mov    %rax,%rdi
    c7e8:	e8 8d f0 ff ff       	call   b87a <_ZN3MPI8DatatypeaSERKP15ompi_datatype_t>
    c7ed:	83 45 f4 01          	addl   $0x1,-0xc(%rbp)
    c7f1:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c7f4:	3b 45 dc             	cmp    -0x24(%rbp),%eax
    c7f7:	7c c2                	jl     c7bb <_ZNK3MPI8Datatype12Get_contentsEiiiPiPlPS0_+0x85>
    c7f9:	48 83 7d f8 00       	cmpq   $0x0,-0x8(%rbp)
    c7fe:	74 0c                	je     c80c <_ZNK3MPI8Datatype12Get_contentsEiiiPiPlPS0_+0xd6>
    c800:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c804:	48 89 c7             	mov    %rax,%rdi
    c807:	e8 b4 dd ff ff       	call   a5c0 <_ZdaPv@plt>
    c80c:	90                   	nop
    c80d:	c9                   	leave  
    c80e:	c3                   	ret    
    c80f:	90                   	nop

000000000000c810 <_ZNK3MPI8Datatype12Get_envelopeERiS1_S1_S1_>:
    c810:	f3 0f 1e fa          	endbr64 
    c814:	55                   	push   %rbp
    c815:	48 89 e5             	mov    %rsp,%rbp
    c818:	48 83 ec 30          	sub    $0x30,%rsp
    c81c:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c820:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c824:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    c828:	48 89 4d e0          	mov    %rcx,-0x20(%rbp)
    c82c:	4c 89 45 d8          	mov    %r8,-0x28(%rbp)
    c830:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c834:	48 8b 40 08          	mov    0x8(%rax),%rax
    c838:	48 8b 7d d8          	mov    -0x28(%rbp),%rdi
    c83c:	48 8b 4d e0          	mov    -0x20(%rbp),%rcx
    c840:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    c844:	48 8b 75 f0          	mov    -0x10(%rbp),%rsi
    c848:	49 89 f8             	mov    %rdi,%r8
    c84b:	48 89 c7             	mov    %rax,%rdi
    c84e:	e8 dd d7 ff ff       	call   a030 <MPI_Type_get_envelope@plt>
    c853:	90                   	nop
    c854:	c9                   	leave  
    c855:	c3                   	ret    

000000000000c856 <_ZNK3MPI8Datatype8Get_nameEPcRi>:
    c856:	f3 0f 1e fa          	endbr64 
    c85a:	55                   	push   %rbp
    c85b:	48 89 e5             	mov    %rsp,%rbp
    c85e:	48 83 ec 20          	sub    $0x20,%rsp
    c862:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c866:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c86a:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    c86e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c872:	48 8b 40 08          	mov    0x8(%rax),%rax
    c876:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    c87a:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    c87e:	48 89 ce             	mov    %rcx,%rsi
    c881:	48 89 c7             	mov    %rax,%rdi
    c884:	e8 17 de ff ff       	call   a6a0 <MPI_Type_get_name@plt>
    c889:	90                   	nop
    c88a:	c9                   	leave  
    c88b:	c3                   	ret    

000000000000c88c <_ZN3MPI8Datatype8Set_attrEiPKv>:
    c88c:	f3 0f 1e fa          	endbr64 
    c890:	55                   	push   %rbp
    c891:	48 89 e5             	mov    %rsp,%rbp
    c894:	48 83 ec 20          	sub    $0x20,%rsp
    c898:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c89c:	89 75 f4             	mov    %esi,-0xc(%rbp)
    c89f:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    c8a3:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c8a7:	48 8b 40 08          	mov    0x8(%rax),%rax
    c8ab:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    c8af:	8b 4d f4             	mov    -0xc(%rbp),%ecx
    c8b2:	89 ce                	mov    %ecx,%esi
    c8b4:	48 89 c7             	mov    %rax,%rdi
    c8b7:	e8 94 e1 ff ff       	call   aa50 <MPI_Type_set_attr@plt>
    c8bc:	90                   	nop
    c8bd:	c9                   	leave  
    c8be:	c3                   	ret    
    c8bf:	90                   	nop

000000000000c8c0 <_ZN3MPI8Datatype8Set_nameEPKc>:
    c8c0:	f3 0f 1e fa          	endbr64 
    c8c4:	55                   	push   %rbp
    c8c5:	48 89 e5             	mov    %rsp,%rbp
    c8c8:	48 83 ec 10          	sub    $0x10,%rsp
    c8cc:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c8d0:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c8d4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c8d8:	48 8b 40 08          	mov    0x8(%rax),%rax
    c8dc:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    c8e0:	48 89 d6             	mov    %rdx,%rsi
    c8e3:	48 89 c7             	mov    %rax,%rdi
    c8e6:	e8 b5 e1 ff ff       	call   aaa0 <MPI_Type_set_name@plt>
    c8eb:	90                   	nop
    c8ec:	c9                   	leave  
    c8ed:	c3                   	ret    

000000000000c8ee <_ZN3MPI14Is_initializedEv>:
    c8ee:	f3 0f 1e fa          	endbr64 
    c8f2:	55                   	push   %rbp
    c8f3:	48 89 e5             	mov    %rsp,%rbp
    c8f6:	48 83 ec 10          	sub    $0x10,%rsp
    c8fa:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c901:	00 00 
    c903:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c907:	31 c0                	xor    %eax,%eax
    c909:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
    c90d:	48 89 c7             	mov    %rax,%rdi
    c910:	e8 4b d7 ff ff       	call   a060 <MPI_Initialized@plt>
    c915:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c918:	85 c0                	test   %eax,%eax
    c91a:	0f 95 c0             	setne  %al
    c91d:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    c921:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    c928:	00 00 
    c92a:	74 05                	je     c931 <_ZN3MPI14Is_initializedEv+0x43>
    c92c:	e8 9f db ff ff       	call   a4d0 <__stack_chk_fail@plt>
    c931:	c9                   	leave  
    c932:	c3                   	ret    
    c933:	90                   	nop

000000000000c934 <_ZN3MPI7Request4WaitERNS_6StatusE>:
    c934:	f3 0f 1e fa          	endbr64 
    c938:	55                   	push   %rbp
    c939:	48 89 e5             	mov    %rsp,%rbp
    c93c:	48 83 ec 10          	sub    $0x10,%rsp
    c940:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c944:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    c948:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    c94c:	48 8d 50 08          	lea    0x8(%rax),%rdx
    c950:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c954:	48 83 c0 08          	add    $0x8,%rax
    c958:	48 89 d6             	mov    %rdx,%rsi
    c95b:	48 89 c7             	mov    %rax,%rdi
    c95e:	e8 cd dc ff ff       	call   a630 <MPI_Wait@plt>
    c963:	90                   	nop
    c964:	c9                   	leave  
    c965:	c3                   	ret    

000000000000c966 <_ZN3MPI7Request4WaitEv>:
    c966:	f3 0f 1e fa          	endbr64 
    c96a:	55                   	push   %rbp
    c96b:	48 89 e5             	mov    %rsp,%rbp
    c96e:	48 83 ec 10          	sub    $0x10,%rsp
    c972:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c976:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c97a:	48 83 c0 08          	add    $0x8,%rax
    c97e:	be 00 00 00 00       	mov    $0x0,%esi
    c983:	48 89 c7             	mov    %rax,%rdi
    c986:	e8 a5 dc ff ff       	call   a630 <MPI_Wait@plt>
    c98b:	90                   	nop
    c98c:	c9                   	leave  
    c98d:	c3                   	ret    

000000000000c98e <_ZN3MPI7Request4FreeEv>:
    c98e:	f3 0f 1e fa          	endbr64 
    c992:	55                   	push   %rbp
    c993:	48 89 e5             	mov    %rsp,%rbp
    c996:	48 83 ec 10          	sub    $0x10,%rsp
    c99a:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    c99e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    c9a2:	48 83 c0 08          	add    $0x8,%rax
    c9a6:	48 89 c7             	mov    %rax,%rdi
    c9a9:	e8 c2 d8 ff ff       	call   a270 <MPI_Request_free@plt>
    c9ae:	90                   	nop
    c9af:	c9                   	leave  
    c9b0:	c3                   	ret    
    c9b1:	90                   	nop

000000000000c9b2 <_ZN3MPI7Request4TestERNS_6StatusE>:
    c9b2:	f3 0f 1e fa          	endbr64 
    c9b6:	55                   	push   %rbp
    c9b7:	48 89 e5             	mov    %rsp,%rbp
    c9ba:	48 83 ec 20          	sub    $0x20,%rsp
    c9be:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    c9c2:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    c9c6:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    c9cd:	00 00 
    c9cf:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    c9d3:	31 c0                	xor    %eax,%eax
    c9d5:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    c9d9:	48 8d 50 08          	lea    0x8(%rax),%rdx
    c9dd:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    c9e1:	48 8d 48 08          	lea    0x8(%rax),%rcx
    c9e5:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
    c9e9:	48 89 c6             	mov    %rax,%rsi
    c9ec:	48 89 cf             	mov    %rcx,%rdi
    c9ef:	e8 bc de ff ff       	call   a8b0 <MPI_Test@plt>
    c9f4:	8b 45 f4             	mov    -0xc(%rbp),%eax
    c9f7:	85 c0                	test   %eax,%eax
    c9f9:	0f 95 c0             	setne  %al
    c9fc:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    ca00:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    ca07:	00 00 
    ca09:	74 05                	je     ca10 <_ZN3MPI7Request4TestERNS_6StatusE+0x5e>
    ca0b:	e8 c0 da ff ff       	call   a4d0 <__stack_chk_fail@plt>
    ca10:	c9                   	leave  
    ca11:	c3                   	ret    

000000000000ca12 <_ZN3MPI7Request4TestEv>:
    ca12:	f3 0f 1e fa          	endbr64 
    ca16:	55                   	push   %rbp
    ca17:	48 89 e5             	mov    %rsp,%rbp
    ca1a:	48 83 ec 20          	sub    $0x20,%rsp
    ca1e:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    ca22:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ca29:	00 00 
    ca2b:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    ca2f:	31 c0                	xor    %eax,%eax
    ca31:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ca35:	48 8d 48 08          	lea    0x8(%rax),%rcx
    ca39:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
    ca3d:	ba 00 00 00 00       	mov    $0x0,%edx
    ca42:	48 89 c6             	mov    %rax,%rsi
    ca45:	48 89 cf             	mov    %rcx,%rdi
    ca48:	e8 63 de ff ff       	call   a8b0 <MPI_Test@plt>
    ca4d:	8b 45 f4             	mov    -0xc(%rbp),%eax
    ca50:	85 c0                	test   %eax,%eax
    ca52:	0f 95 c0             	setne  %al
    ca55:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    ca59:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    ca60:	00 00 
    ca62:	74 05                	je     ca69 <_ZN3MPI7Request4TestEv+0x57>
    ca64:	e8 67 da ff ff       	call   a4d0 <__stack_chk_fail@plt>
    ca69:	c9                   	leave  
    ca6a:	c3                   	ret    
    ca6b:	90                   	nop

000000000000ca6c <_ZNK3MPI7Request6CancelEv>:
    ca6c:	f3 0f 1e fa          	endbr64 
    ca70:	55                   	push   %rbp
    ca71:	48 89 e5             	mov    %rsp,%rbp
    ca74:	48 83 ec 10          	sub    $0x10,%rsp
    ca78:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    ca7c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ca80:	48 83 c0 08          	add    $0x8,%rax
    ca84:	48 89 c7             	mov    %rax,%rdi
    ca87:	e8 c4 d5 ff ff       	call   a050 <MPI_Cancel@plt>
    ca8c:	90                   	nop
    ca8d:	c9                   	leave  
    ca8e:	c3                   	ret    
    ca8f:	90                   	nop

000000000000ca90 <_ZN3MPI8Prequest5StartEv>:
    ca90:	f3 0f 1e fa          	endbr64 
    ca94:	55                   	push   %rbp
    ca95:	48 89 e5             	mov    %rsp,%rbp
    ca98:	48 83 ec 10          	sub    $0x10,%rsp
    ca9c:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    caa0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    caa4:	48 83 c0 08          	add    $0x8,%rax
    caa8:	48 89 c7             	mov    %rax,%rdi
    caab:	e8 40 de ff ff       	call   a8f0 <MPI_Start@plt>
    cab0:	90                   	nop
    cab1:	c9                   	leave  
    cab2:	c3                   	ret    
    cab3:	90                   	nop

000000000000cab4 <_ZNK3MPI7Request10Get_statusERNS_6StatusE>:
    cab4:	f3 0f 1e fa          	endbr64 
    cab8:	55                   	push   %rbp
    cab9:	48 89 e5             	mov    %rsp,%rbp
    cabc:	48 83 ec 40          	sub    $0x40,%rsp
    cac0:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
    cac4:	48 89 75 c0          	mov    %rsi,-0x40(%rbp)
    cac8:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    cacf:	00 00 
    cad1:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    cad5:	31 c0                	xor    %eax,%eax
    cad7:	c7 45 dc 00 00 00 00 	movl   $0x0,-0x24(%rbp)
    cade:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    cae2:	48 8b 40 08          	mov    0x8(%rax),%rax
    cae6:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    caea:	48 8d 4d dc          	lea    -0x24(%rbp),%rcx
    caee:	48 89 ce             	mov    %rcx,%rsi
    caf1:	48 89 c7             	mov    %rax,%rdi
    caf4:	e8 27 d4 ff ff       	call   9f20 <MPI_Request_get_status@plt>
    caf9:	8b 45 dc             	mov    -0x24(%rbp),%eax
    cafc:	85 c0                	test   %eax,%eax
    cafe:	74 13                	je     cb13 <_ZNK3MPI7Request10Get_statusERNS_6StatusE+0x5f>
    cb00:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    cb04:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    cb08:	48 89 d6             	mov    %rdx,%rsi
    cb0b:	48 89 c7             	mov    %rax,%rdi
    cb0e:	e8 f1 ed ff ff       	call   b904 <_ZN3MPI6StatusaSERK20ompi_status_public_t>
    cb13:	8b 45 dc             	mov    -0x24(%rbp),%eax
    cb16:	85 c0                	test   %eax,%eax
    cb18:	0f 95 c0             	setne  %al
    cb1b:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    cb1f:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    cb26:	00 00 
    cb28:	74 05                	je     cb2f <_ZNK3MPI7Request10Get_statusERNS_6StatusE+0x7b>
    cb2a:	e8 a1 d9 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    cb2f:	c9                   	leave  
    cb30:	c3                   	ret    
    cb31:	90                   	nop

000000000000cb32 <_ZNK3MPI7Request10Get_statusEv>:
    cb32:	f3 0f 1e fa          	endbr64 
    cb36:	55                   	push   %rbp
    cb37:	48 89 e5             	mov    %rsp,%rbp
    cb3a:	48 83 ec 20          	sub    $0x20,%rsp
    cb3e:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cb42:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    cb49:	00 00 
    cb4b:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    cb4f:	31 c0                	xor    %eax,%eax
    cb51:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cb55:	48 8b 40 08          	mov    0x8(%rax),%rax
    cb59:	48 8d 4d f4          	lea    -0xc(%rbp),%rcx
    cb5d:	ba 00 00 00 00       	mov    $0x0,%edx
    cb62:	48 89 ce             	mov    %rcx,%rsi
    cb65:	48 89 c7             	mov    %rax,%rdi
    cb68:	e8 b3 d3 ff ff       	call   9f20 <MPI_Request_get_status@plt>
    cb6d:	8b 45 f4             	mov    -0xc(%rbp),%eax
    cb70:	85 c0                	test   %eax,%eax
    cb72:	0f 95 c0             	setne  %al
    cb75:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    cb79:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    cb80:	00 00 
    cb82:	74 05                	je     cb89 <_ZNK3MPI7Request10Get_statusEv+0x57>
    cb84:	e8 47 d9 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    cb89:	c9                   	leave  
    cb8a:	c3                   	ret    
    cb8b:	90                   	nop

000000000000cb8c <_ZN3MPI8Grequest8CompleteEv>:
    cb8c:	f3 0f 1e fa          	endbr64 
    cb90:	55                   	push   %rbp
    cb91:	48 89 e5             	mov    %rsp,%rbp
    cb94:	48 83 ec 10          	sub    $0x10,%rsp
    cb98:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    cb9c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    cba0:	48 8b 40 08          	mov    0x8(%rax),%rax
    cba4:	48 89 c7             	mov    %rax,%rdi
    cba7:	e8 e4 d6 ff ff       	call   a290 <MPI_Grequest_complete@plt>
    cbac:	90                   	nop
    cbad:	c9                   	leave  
    cbae:	c3                   	ret    
    cbaf:	90                   	nop

000000000000cbb0 <_ZNK3MPI4Comm4SendEPKviRKNS_8DatatypeEii>:
    cbb0:	f3 0f 1e fa          	endbr64 
    cbb4:	55                   	push   %rbp
    cbb5:	48 89 e5             	mov    %rsp,%rbp
    cbb8:	53                   	push   %rbx
    cbb9:	48 83 ec 38          	sub    $0x38,%rsp
    cbbd:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cbc1:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    cbc5:	89 55 dc             	mov    %edx,-0x24(%rbp)
    cbc8:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    cbcc:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    cbd0:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
    cbd4:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cbd8:	48 8b 58 08          	mov    0x8(%rax),%rbx
    cbdc:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cbe0:	48 89 c7             	mov    %rax,%rdi
    cbe3:	e8 b8 ec ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cbe8:	48 89 c7             	mov    %rax,%rdi
    cbeb:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    cbee:	8b 55 d8             	mov    -0x28(%rbp),%edx
    cbf1:	8b 75 dc             	mov    -0x24(%rbp),%esi
    cbf4:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    cbf8:	49 89 d9             	mov    %rbx,%r9
    cbfb:	41 89 c8             	mov    %ecx,%r8d
    cbfe:	89 d1                	mov    %edx,%ecx
    cc00:	48 89 fa             	mov    %rdi,%rdx
    cc03:	48 89 c7             	mov    %rax,%rdi
    cc06:	e8 75 d3 ff ff       	call   9f80 <MPI_Send@plt>
    cc0b:	90                   	nop
    cc0c:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    cc10:	c9                   	leave  
    cc11:	c3                   	ret    

000000000000cc12 <_ZNK3MPI4Comm4RecvEPviRKNS_8DatatypeEiiRNS_6StatusE>:
    cc12:	f3 0f 1e fa          	endbr64 
    cc16:	55                   	push   %rbp
    cc17:	48 89 e5             	mov    %rsp,%rbp
    cc1a:	41 54                	push   %r12
    cc1c:	53                   	push   %rbx
    cc1d:	48 83 ec 30          	sub    $0x30,%rsp
    cc21:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cc25:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    cc29:	89 55 dc             	mov    %edx,-0x24(%rbp)
    cc2c:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    cc30:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    cc34:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
    cc38:	48 8b 45 10          	mov    0x10(%rbp),%rax
    cc3c:	48 8d 58 08          	lea    0x8(%rax),%rbx
    cc40:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cc44:	4c 8b 60 08          	mov    0x8(%rax),%r12
    cc48:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cc4c:	48 89 c7             	mov    %rax,%rdi
    cc4f:	e8 4c ec ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cc54:	48 89 c7             	mov    %rax,%rdi
    cc57:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    cc5a:	8b 55 d8             	mov    -0x28(%rbp),%edx
    cc5d:	8b 75 dc             	mov    -0x24(%rbp),%esi
    cc60:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    cc64:	48 83 ec 08          	sub    $0x8,%rsp
    cc68:	53                   	push   %rbx
    cc69:	4d 89 e1             	mov    %r12,%r9
    cc6c:	41 89 c8             	mov    %ecx,%r8d
    cc6f:	89 d1                	mov    %edx,%ecx
    cc71:	48 89 fa             	mov    %rdi,%rdx
    cc74:	48 89 c7             	mov    %rax,%rdi
    cc77:	e8 54 d6 ff ff       	call   a2d0 <MPI_Recv@plt>
    cc7c:	48 83 c4 10          	add    $0x10,%rsp
    cc80:	90                   	nop
    cc81:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    cc85:	5b                   	pop    %rbx
    cc86:	41 5c                	pop    %r12
    cc88:	5d                   	pop    %rbp
    cc89:	c3                   	ret    

000000000000cc8a <_ZNK3MPI4Comm4RecvEPviRKNS_8DatatypeEii>:
    cc8a:	f3 0f 1e fa          	endbr64 
    cc8e:	55                   	push   %rbp
    cc8f:	48 89 e5             	mov    %rsp,%rbp
    cc92:	53                   	push   %rbx
    cc93:	48 83 ec 38          	sub    $0x38,%rsp
    cc97:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cc9b:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    cc9f:	89 55 dc             	mov    %edx,-0x24(%rbp)
    cca2:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    cca6:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    ccaa:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
    ccae:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ccb2:	48 8b 58 08          	mov    0x8(%rax),%rbx
    ccb6:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ccba:	48 89 c7             	mov    %rax,%rdi
    ccbd:	e8 de eb ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    ccc2:	48 89 c7             	mov    %rax,%rdi
    ccc5:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    ccc8:	8b 55 d8             	mov    -0x28(%rbp),%edx
    cccb:	8b 75 dc             	mov    -0x24(%rbp),%esi
    ccce:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    ccd2:	48 83 ec 08          	sub    $0x8,%rsp
    ccd6:	6a 00                	push   $0x0
    ccd8:	49 89 d9             	mov    %rbx,%r9
    ccdb:	41 89 c8             	mov    %ecx,%r8d
    ccde:	89 d1                	mov    %edx,%ecx
    cce0:	48 89 fa             	mov    %rdi,%rdx
    cce3:	48 89 c7             	mov    %rax,%rdi
    cce6:	e8 e5 d5 ff ff       	call   a2d0 <MPI_Recv@plt>
    cceb:	48 83 c4 10          	add    $0x10,%rsp
    ccef:	90                   	nop
    ccf0:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    ccf4:	c9                   	leave  
    ccf5:	c3                   	ret    

000000000000ccf6 <_ZNK3MPI4Comm5BsendEPKviRKNS_8DatatypeEii>:
    ccf6:	f3 0f 1e fa          	endbr64 
    ccfa:	55                   	push   %rbp
    ccfb:	48 89 e5             	mov    %rsp,%rbp
    ccfe:	53                   	push   %rbx
    ccff:	48 83 ec 38          	sub    $0x38,%rsp
    cd03:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cd07:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    cd0b:	89 55 dc             	mov    %edx,-0x24(%rbp)
    cd0e:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    cd12:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    cd16:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
    cd1a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cd1e:	48 8b 58 08          	mov    0x8(%rax),%rbx
    cd22:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cd26:	48 89 c7             	mov    %rax,%rdi
    cd29:	e8 72 eb ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cd2e:	48 89 c7             	mov    %rax,%rdi
    cd31:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    cd34:	8b 55 d8             	mov    -0x28(%rbp),%edx
    cd37:	8b 75 dc             	mov    -0x24(%rbp),%esi
    cd3a:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    cd3e:	49 89 d9             	mov    %rbx,%r9
    cd41:	41 89 c8             	mov    %ecx,%r8d
    cd44:	89 d1                	mov    %edx,%ecx
    cd46:	48 89 fa             	mov    %rdi,%rdx
    cd49:	48 89 c7             	mov    %rax,%rdi
    cd4c:	e8 5f d7 ff ff       	call   a4b0 <MPI_Bsend@plt>
    cd51:	90                   	nop
    cd52:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    cd56:	c9                   	leave  
    cd57:	c3                   	ret    

000000000000cd58 <_ZNK3MPI4Comm5SsendEPKviRKNS_8DatatypeEii>:
    cd58:	f3 0f 1e fa          	endbr64 
    cd5c:	55                   	push   %rbp
    cd5d:	48 89 e5             	mov    %rsp,%rbp
    cd60:	53                   	push   %rbx
    cd61:	48 83 ec 38          	sub    $0x38,%rsp
    cd65:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cd69:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    cd6d:	89 55 dc             	mov    %edx,-0x24(%rbp)
    cd70:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    cd74:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    cd78:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
    cd7c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cd80:	48 8b 58 08          	mov    0x8(%rax),%rbx
    cd84:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cd88:	48 89 c7             	mov    %rax,%rdi
    cd8b:	e8 10 eb ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cd90:	48 89 c7             	mov    %rax,%rdi
    cd93:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    cd96:	8b 55 d8             	mov    -0x28(%rbp),%edx
    cd99:	8b 75 dc             	mov    -0x24(%rbp),%esi
    cd9c:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    cda0:	49 89 d9             	mov    %rbx,%r9
    cda3:	41 89 c8             	mov    %ecx,%r8d
    cda6:	89 d1                	mov    %edx,%ecx
    cda8:	48 89 fa             	mov    %rdi,%rdx
    cdab:	48 89 c7             	mov    %rax,%rdi
    cdae:	e8 bd d0 ff ff       	call   9e70 <MPI_Ssend@plt>
    cdb3:	90                   	nop
    cdb4:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    cdb8:	c9                   	leave  
    cdb9:	c3                   	ret    

000000000000cdba <_ZNK3MPI4Comm5RsendEPKviRKNS_8DatatypeEii>:
    cdba:	f3 0f 1e fa          	endbr64 
    cdbe:	55                   	push   %rbp
    cdbf:	48 89 e5             	mov    %rsp,%rbp
    cdc2:	53                   	push   %rbx
    cdc3:	48 83 ec 38          	sub    $0x38,%rsp
    cdc7:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    cdcb:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    cdcf:	89 55 dc             	mov    %edx,-0x24(%rbp)
    cdd2:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    cdd6:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    cdda:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
    cdde:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cde2:	48 8b 58 08          	mov    0x8(%rax),%rbx
    cde6:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cdea:	48 89 c7             	mov    %rax,%rdi
    cded:	e8 ae ea ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cdf2:	48 89 c7             	mov    %rax,%rdi
    cdf5:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    cdf8:	8b 55 d8             	mov    -0x28(%rbp),%edx
    cdfb:	8b 75 dc             	mov    -0x24(%rbp),%esi
    cdfe:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    ce02:	49 89 d9             	mov    %rbx,%r9
    ce05:	41 89 c8             	mov    %ecx,%r8d
    ce08:	89 d1                	mov    %edx,%ecx
    ce0a:	48 89 fa             	mov    %rdi,%rdx
    ce0d:	48 89 c7             	mov    %rax,%rdi
    ce10:	e8 ab d4 ff ff       	call   a2c0 <MPI_Rsend@plt>
    ce15:	90                   	nop
    ce16:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    ce1a:	c9                   	leave  
    ce1b:	c3                   	ret    

000000000000ce1c <_ZNK3MPI4Comm5IsendEPKviRKNS_8DatatypeEii>:
    ce1c:	f3 0f 1e fa          	endbr64 
    ce20:	55                   	push   %rbp
    ce21:	48 89 e5             	mov    %rsp,%rbp
    ce24:	53                   	push   %rbx
    ce25:	48 83 ec 48          	sub    $0x48,%rsp
    ce29:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    ce2d:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    ce31:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    ce35:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    ce38:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    ce3c:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    ce40:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ce47:	00 00 
    ce49:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    ce4d:	31 c0                	xor    %eax,%eax
    ce4f:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ce53:	48 8b 58 08          	mov    0x8(%rax),%rbx
    ce57:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    ce5b:	48 89 c7             	mov    %rax,%rdi
    ce5e:	e8 3d ea ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    ce63:	48 89 c7             	mov    %rax,%rdi
    ce66:	8b 55 c0             	mov    -0x40(%rbp),%edx
    ce69:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    ce6c:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    ce70:	48 83 ec 08          	sub    $0x8,%rsp
    ce74:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    ce78:	51                   	push   %rcx
    ce79:	49 89 d9             	mov    %rbx,%r9
    ce7c:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    ce80:	89 d1                	mov    %edx,%ecx
    ce82:	48 89 fa             	mov    %rdi,%rdx
    ce85:	48 89 c7             	mov    %rax,%rdi
    ce88:	e8 63 d1 ff ff       	call   9ff0 <MPI_Isend@plt>
    ce8d:	48 83 c4 10          	add    $0x10,%rsp
    ce91:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    ce95:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ce99:	48 89 d6             	mov    %rdx,%rsi
    ce9c:	48 89 c7             	mov    %rax,%rdi
    ce9f:	e8 e4 ea ff ff       	call   b988 <_ZN3MPI7RequestC1EP14ompi_request_t>
    cea4:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cea8:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    ceaf:	00 00 
    ceb1:	74 05                	je     ceb8 <_ZNK3MPI4Comm5IsendEPKviRKNS_8DatatypeEii+0x9c>
    ceb3:	e8 18 d6 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    ceb8:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    cebc:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    cec0:	c9                   	leave  
    cec1:	c3                   	ret    

000000000000cec2 <_ZNK3MPI4Comm6IbsendEPKviRKNS_8DatatypeEii>:
    cec2:	f3 0f 1e fa          	endbr64 
    cec6:	55                   	push   %rbp
    cec7:	48 89 e5             	mov    %rsp,%rbp
    ceca:	53                   	push   %rbx
    cecb:	48 83 ec 48          	sub    $0x48,%rsp
    cecf:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    ced3:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    ced7:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    cedb:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    cede:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    cee2:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    cee6:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ceed:	00 00 
    ceef:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    cef3:	31 c0                	xor    %eax,%eax
    cef5:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cef9:	48 8b 58 08          	mov    0x8(%rax),%rbx
    cefd:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    cf01:	48 89 c7             	mov    %rax,%rdi
    cf04:	e8 97 e9 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cf09:	48 89 c7             	mov    %rax,%rdi
    cf0c:	8b 55 c0             	mov    -0x40(%rbp),%edx
    cf0f:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    cf12:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    cf16:	48 83 ec 08          	sub    $0x8,%rsp
    cf1a:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    cf1e:	51                   	push   %rcx
    cf1f:	49 89 d9             	mov    %rbx,%r9
    cf22:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    cf26:	89 d1                	mov    %edx,%ecx
    cf28:	48 89 fa             	mov    %rdi,%rdx
    cf2b:	48 89 c7             	mov    %rax,%rdi
    cf2e:	e8 ad d8 ff ff       	call   a7e0 <MPI_Ibsend@plt>
    cf33:	48 83 c4 10          	add    $0x10,%rsp
    cf37:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    cf3b:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    cf3f:	48 89 d6             	mov    %rdx,%rsi
    cf42:	48 89 c7             	mov    %rax,%rdi
    cf45:	e8 3e ea ff ff       	call   b988 <_ZN3MPI7RequestC1EP14ompi_request_t>
    cf4a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cf4e:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    cf55:	00 00 
    cf57:	74 05                	je     cf5e <_ZNK3MPI4Comm6IbsendEPKviRKNS_8DatatypeEii+0x9c>
    cf59:	e8 72 d5 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    cf5e:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    cf62:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    cf66:	c9                   	leave  
    cf67:	c3                   	ret    

000000000000cf68 <_ZNK3MPI4Comm6IssendEPKviRKNS_8DatatypeEii>:
    cf68:	f3 0f 1e fa          	endbr64 
    cf6c:	55                   	push   %rbp
    cf6d:	48 89 e5             	mov    %rsp,%rbp
    cf70:	53                   	push   %rbx
    cf71:	48 83 ec 48          	sub    $0x48,%rsp
    cf75:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    cf79:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    cf7d:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    cf81:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    cf84:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    cf88:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    cf8c:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    cf93:	00 00 
    cf95:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    cf99:	31 c0                	xor    %eax,%eax
    cf9b:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    cf9f:	48 8b 58 08          	mov    0x8(%rax),%rbx
    cfa3:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    cfa7:	48 89 c7             	mov    %rax,%rdi
    cfaa:	e8 f1 e8 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    cfaf:	48 89 c7             	mov    %rax,%rdi
    cfb2:	8b 55 c0             	mov    -0x40(%rbp),%edx
    cfb5:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    cfb8:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    cfbc:	48 83 ec 08          	sub    $0x8,%rsp
    cfc0:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    cfc4:	51                   	push   %rcx
    cfc5:	49 89 d9             	mov    %rbx,%r9
    cfc8:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    cfcc:	89 d1                	mov    %edx,%ecx
    cfce:	48 89 fa             	mov    %rdi,%rdx
    cfd1:	48 89 c7             	mov    %rax,%rdi
    cfd4:	e8 a7 d0 ff ff       	call   a080 <MPI_Issend@plt>
    cfd9:	48 83 c4 10          	add    $0x10,%rsp
    cfdd:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    cfe1:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    cfe5:	48 89 d6             	mov    %rdx,%rsi
    cfe8:	48 89 c7             	mov    %rax,%rdi
    cfeb:	e8 98 e9 ff ff       	call   b988 <_ZN3MPI7RequestC1EP14ompi_request_t>
    cff0:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    cff4:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    cffb:	00 00 
    cffd:	74 05                	je     d004 <_ZNK3MPI4Comm6IssendEPKviRKNS_8DatatypeEii+0x9c>
    cfff:	e8 cc d4 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d004:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d008:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d00c:	c9                   	leave  
    d00d:	c3                   	ret    

000000000000d00e <_ZNK3MPI4Comm6IrsendEPKviRKNS_8DatatypeEii>:
    d00e:	f3 0f 1e fa          	endbr64 
    d012:	55                   	push   %rbp
    d013:	48 89 e5             	mov    %rsp,%rbp
    d016:	53                   	push   %rbx
    d017:	48 83 ec 48          	sub    $0x48,%rsp
    d01b:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d01f:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d023:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d027:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d02a:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d02e:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d032:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d039:	00 00 
    d03b:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d03f:	31 c0                	xor    %eax,%eax
    d041:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d045:	48 8b 58 08          	mov    0x8(%rax),%rbx
    d049:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d04d:	48 89 c7             	mov    %rax,%rdi
    d050:	e8 4b e8 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d055:	48 89 c7             	mov    %rax,%rdi
    d058:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d05b:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d05e:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d062:	48 83 ec 08          	sub    $0x8,%rsp
    d066:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d06a:	51                   	push   %rcx
    d06b:	49 89 d9             	mov    %rbx,%r9
    d06e:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d072:	89 d1                	mov    %edx,%ecx
    d074:	48 89 fa             	mov    %rdi,%rdx
    d077:	48 89 c7             	mov    %rax,%rdi
    d07a:	e8 c1 cd ff ff       	call   9e40 <MPI_Irsend@plt>
    d07f:	48 83 c4 10          	add    $0x10,%rsp
    d083:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    d087:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d08b:	48 89 d6             	mov    %rdx,%rsi
    d08e:	48 89 c7             	mov    %rax,%rdi
    d091:	e8 f2 e8 ff ff       	call   b988 <_ZN3MPI7RequestC1EP14ompi_request_t>
    d096:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d09a:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d0a1:	00 00 
    d0a3:	74 05                	je     d0aa <_ZNK3MPI4Comm6IrsendEPKviRKNS_8DatatypeEii+0x9c>
    d0a5:	e8 26 d4 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d0aa:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d0ae:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d0b2:	c9                   	leave  
    d0b3:	c3                   	ret    

000000000000d0b4 <_ZNK3MPI4Comm5IrecvEPviRKNS_8DatatypeEii>:
    d0b4:	f3 0f 1e fa          	endbr64 
    d0b8:	55                   	push   %rbp
    d0b9:	48 89 e5             	mov    %rsp,%rbp
    d0bc:	53                   	push   %rbx
    d0bd:	48 83 ec 48          	sub    $0x48,%rsp
    d0c1:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d0c5:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d0c9:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d0cd:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d0d0:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d0d4:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d0d8:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d0df:	00 00 
    d0e1:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d0e5:	31 c0                	xor    %eax,%eax
    d0e7:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d0eb:	48 8b 58 08          	mov    0x8(%rax),%rbx
    d0ef:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d0f3:	48 89 c7             	mov    %rax,%rdi
    d0f6:	e8 a5 e7 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d0fb:	48 89 c7             	mov    %rax,%rdi
    d0fe:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d101:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d104:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d108:	48 83 ec 08          	sub    $0x8,%rsp
    d10c:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d110:	51                   	push   %rcx
    d111:	49 89 d9             	mov    %rbx,%r9
    d114:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d118:	89 d1                	mov    %edx,%ecx
    d11a:	48 89 fa             	mov    %rdi,%rdx
    d11d:	48 89 c7             	mov    %rax,%rdi
    d120:	e8 4b d5 ff ff       	call   a670 <MPI_Irecv@plt>
    d125:	48 83 c4 10          	add    $0x10,%rsp
    d129:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    d12d:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d131:	48 89 d6             	mov    %rdx,%rsi
    d134:	48 89 c7             	mov    %rax,%rdi
    d137:	e8 4c e8 ff ff       	call   b988 <_ZN3MPI7RequestC1EP14ompi_request_t>
    d13c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d140:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d147:	00 00 
    d149:	74 05                	je     d150 <_ZNK3MPI4Comm5IrecvEPviRKNS_8DatatypeEii+0x9c>
    d14b:	e8 80 d3 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d150:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d154:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d158:	c9                   	leave  
    d159:	c3                   	ret    

000000000000d15a <_ZNK3MPI4Comm6IprobeEiiRNS_6StatusE>:
    d15a:	f3 0f 1e fa          	endbr64 
    d15e:	55                   	push   %rbp
    d15f:	48 89 e5             	mov    %rsp,%rbp
    d162:	48 83 ec 30          	sub    $0x30,%rsp
    d166:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d16a:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    d16d:	89 55 e0             	mov    %edx,-0x20(%rbp)
    d170:	48 89 4d d8          	mov    %rcx,-0x28(%rbp)
    d174:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d17b:	00 00 
    d17d:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    d181:	31 c0                	xor    %eax,%eax
    d183:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d187:	48 8d 78 08          	lea    0x8(%rax),%rdi
    d18b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d18f:	48 8b 50 08          	mov    0x8(%rax),%rdx
    d193:	48 8d 4d f4          	lea    -0xc(%rbp),%rcx
    d197:	8b 75 e0             	mov    -0x20(%rbp),%esi
    d19a:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    d19d:	49 89 f8             	mov    %rdi,%r8
    d1a0:	89 c7                	mov    %eax,%edi
    d1a2:	e8 69 ce ff ff       	call   a010 <MPI_Iprobe@plt>
    d1a7:	8b 45 f4             	mov    -0xc(%rbp),%eax
    d1aa:	85 c0                	test   %eax,%eax
    d1ac:	0f 95 c0             	setne  %al
    d1af:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    d1b3:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    d1ba:	00 00 
    d1bc:	74 05                	je     d1c3 <_ZNK3MPI4Comm6IprobeEiiRNS_6StatusE+0x69>
    d1be:	e8 0d d3 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d1c3:	c9                   	leave  
    d1c4:	c3                   	ret    
    d1c5:	90                   	nop

000000000000d1c6 <_ZNK3MPI4Comm6IprobeEii>:
    d1c6:	f3 0f 1e fa          	endbr64 
    d1ca:	55                   	push   %rbp
    d1cb:	48 89 e5             	mov    %rsp,%rbp
    d1ce:	48 83 ec 20          	sub    $0x20,%rsp
    d1d2:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d1d6:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    d1d9:	89 55 e0             	mov    %edx,-0x20(%rbp)
    d1dc:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d1e3:	00 00 
    d1e5:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    d1e9:	31 c0                	xor    %eax,%eax
    d1eb:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d1ef:	48 8b 50 08          	mov    0x8(%rax),%rdx
    d1f3:	48 8d 4d f4          	lea    -0xc(%rbp),%rcx
    d1f7:	8b 75 e0             	mov    -0x20(%rbp),%esi
    d1fa:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    d1fd:	41 b8 00 00 00 00    	mov    $0x0,%r8d
    d203:	89 c7                	mov    %eax,%edi
    d205:	e8 06 ce ff ff       	call   a010 <MPI_Iprobe@plt>
    d20a:	8b 45 f4             	mov    -0xc(%rbp),%eax
    d20d:	85 c0                	test   %eax,%eax
    d20f:	0f 95 c0             	setne  %al
    d212:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    d216:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    d21d:	00 00 
    d21f:	74 05                	je     d226 <_ZNK3MPI4Comm6IprobeEii+0x60>
    d221:	e8 aa d2 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d226:	c9                   	leave  
    d227:	c3                   	ret    

000000000000d228 <_ZNK3MPI4Comm5ProbeEiiRNS_6StatusE>:
    d228:	f3 0f 1e fa          	endbr64 
    d22c:	55                   	push   %rbp
    d22d:	48 89 e5             	mov    %rsp,%rbp
    d230:	48 83 ec 20          	sub    $0x20,%rsp
    d234:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    d238:	89 75 f4             	mov    %esi,-0xc(%rbp)
    d23b:	89 55 f0             	mov    %edx,-0x10(%rbp)
    d23e:	48 89 4d e8          	mov    %rcx,-0x18(%rbp)
    d242:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d246:	48 8d 48 08          	lea    0x8(%rax),%rcx
    d24a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    d24e:	48 8b 50 08          	mov    0x8(%rax),%rdx
    d252:	8b 75 f0             	mov    -0x10(%rbp),%esi
    d255:	8b 45 f4             	mov    -0xc(%rbp),%eax
    d258:	89 c7                	mov    %eax,%edi
    d25a:	e8 01 cd ff ff       	call   9f60 <MPI_Probe@plt>
    d25f:	90                   	nop
    d260:	c9                   	leave  
    d261:	c3                   	ret    

000000000000d262 <_ZNK3MPI4Comm5ProbeEii>:
    d262:	f3 0f 1e fa          	endbr64 
    d266:	55                   	push   %rbp
    d267:	48 89 e5             	mov    %rsp,%rbp
    d26a:	48 83 ec 10          	sub    $0x10,%rsp
    d26e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    d272:	89 75 f4             	mov    %esi,-0xc(%rbp)
    d275:	89 55 f0             	mov    %edx,-0x10(%rbp)
    d278:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    d27c:	48 8b 50 08          	mov    0x8(%rax),%rdx
    d280:	8b 75 f0             	mov    -0x10(%rbp),%esi
    d283:	8b 45 f4             	mov    -0xc(%rbp),%eax
    d286:	b9 00 00 00 00       	mov    $0x0,%ecx
    d28b:	89 c7                	mov    %eax,%edi
    d28d:	e8 ce cc ff ff       	call   9f60 <MPI_Probe@plt>
    d292:	90                   	nop
    d293:	c9                   	leave  
    d294:	c3                   	ret    
    d295:	90                   	nop

000000000000d296 <_ZNK3MPI4Comm9Send_initEPKviRKNS_8DatatypeEii>:
    d296:	f3 0f 1e fa          	endbr64 
    d29a:	55                   	push   %rbp
    d29b:	48 89 e5             	mov    %rsp,%rbp
    d29e:	53                   	push   %rbx
    d29f:	48 83 ec 48          	sub    $0x48,%rsp
    d2a3:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d2a7:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d2ab:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d2af:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d2b2:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d2b6:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d2ba:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d2c1:	00 00 
    d2c3:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d2c7:	31 c0                	xor    %eax,%eax
    d2c9:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d2cd:	48 8b 58 08          	mov    0x8(%rax),%rbx
    d2d1:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d2d5:	48 89 c7             	mov    %rax,%rdi
    d2d8:	e8 c3 e5 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d2dd:	48 89 c7             	mov    %rax,%rdi
    d2e0:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d2e3:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d2e6:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d2ea:	48 83 ec 08          	sub    $0x8,%rsp
    d2ee:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d2f2:	51                   	push   %rcx
    d2f3:	49 89 d9             	mov    %rbx,%r9
    d2f6:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d2fa:	89 d1                	mov    %edx,%ecx
    d2fc:	48 89 fa             	mov    %rdi,%rdx
    d2ff:	48 89 c7             	mov    %rax,%rdi
    d302:	e8 89 d1 ff ff       	call   a490 <MPI_Send_init@plt>
    d307:	48 83 c4 10          	add    $0x10,%rsp
    d30b:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    d30f:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d313:	48 89 d6             	mov    %rdx,%rsi
    d316:	48 89 c7             	mov    %rax,%rdi
    d319:	e8 98 e6 ff ff       	call   b9b6 <_ZN3MPI8PrequestC1ERKP14ompi_request_t>
    d31e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d322:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d329:	00 00 
    d32b:	74 05                	je     d332 <_ZNK3MPI4Comm9Send_initEPKviRKNS_8DatatypeEii+0x9c>
    d32d:	e8 9e d1 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d332:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d336:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d33a:	c9                   	leave  
    d33b:	c3                   	ret    

000000000000d33c <_ZNK3MPI4Comm10Bsend_initEPKviRKNS_8DatatypeEii>:
    d33c:	f3 0f 1e fa          	endbr64 
    d340:	55                   	push   %rbp
    d341:	48 89 e5             	mov    %rsp,%rbp
    d344:	53                   	push   %rbx
    d345:	48 83 ec 48          	sub    $0x48,%rsp
    d349:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d34d:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d351:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d355:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d358:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d35c:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d360:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d367:	00 00 
    d369:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d36d:	31 c0                	xor    %eax,%eax
    d36f:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d373:	48 8b 58 08          	mov    0x8(%rax),%rbx
    d377:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d37b:	48 89 c7             	mov    %rax,%rdi
    d37e:	e8 1d e5 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d383:	48 89 c7             	mov    %rax,%rdi
    d386:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d389:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d38c:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d390:	48 83 ec 08          	sub    $0x8,%rsp
    d394:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d398:	51                   	push   %rcx
    d399:	49 89 d9             	mov    %rbx,%r9
    d39c:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d3a0:	89 d1                	mov    %edx,%ecx
    d3a2:	48 89 fa             	mov    %rdi,%rdx
    d3a5:	48 89 c7             	mov    %rax,%rdi
    d3a8:	e8 73 ce ff ff       	call   a220 <MPI_Bsend_init@plt>
    d3ad:	48 83 c4 10          	add    $0x10,%rsp
    d3b1:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    d3b5:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d3b9:	48 89 d6             	mov    %rdx,%rsi
    d3bc:	48 89 c7             	mov    %rax,%rdi
    d3bf:	e8 f2 e5 ff ff       	call   b9b6 <_ZN3MPI8PrequestC1ERKP14ompi_request_t>
    d3c4:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d3c8:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d3cf:	00 00 
    d3d1:	74 05                	je     d3d8 <_ZNK3MPI4Comm10Bsend_initEPKviRKNS_8DatatypeEii+0x9c>
    d3d3:	e8 f8 d0 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d3d8:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d3dc:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d3e0:	c9                   	leave  
    d3e1:	c3                   	ret    

000000000000d3e2 <_ZNK3MPI4Comm10Ssend_initEPKviRKNS_8DatatypeEii>:
    d3e2:	f3 0f 1e fa          	endbr64 
    d3e6:	55                   	push   %rbp
    d3e7:	48 89 e5             	mov    %rsp,%rbp
    d3ea:	53                   	push   %rbx
    d3eb:	48 83 ec 48          	sub    $0x48,%rsp
    d3ef:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d3f3:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d3f7:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d3fb:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d3fe:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d402:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d406:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d40d:	00 00 
    d40f:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d413:	31 c0                	xor    %eax,%eax
    d415:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d419:	48 8b 58 08          	mov    0x8(%rax),%rbx
    d41d:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d421:	48 89 c7             	mov    %rax,%rdi
    d424:	e8 77 e4 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d429:	48 89 c7             	mov    %rax,%rdi
    d42c:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d42f:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d432:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d436:	48 83 ec 08          	sub    $0x8,%rsp
    d43a:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d43e:	51                   	push   %rcx
    d43f:	49 89 d9             	mov    %rbx,%r9
    d442:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d446:	89 d1                	mov    %edx,%ecx
    d448:	48 89 fa             	mov    %rdi,%rdx
    d44b:	48 89 c7             	mov    %rax,%rdi
    d44e:	e8 2d cd ff ff       	call   a180 <MPI_Ssend_init@plt>
    d453:	48 83 c4 10          	add    $0x10,%rsp
    d457:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    d45b:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d45f:	48 89 d6             	mov    %rdx,%rsi
    d462:	48 89 c7             	mov    %rax,%rdi
    d465:	e8 4c e5 ff ff       	call   b9b6 <_ZN3MPI8PrequestC1ERKP14ompi_request_t>
    d46a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d46e:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d475:	00 00 
    d477:	74 05                	je     d47e <_ZNK3MPI4Comm10Ssend_initEPKviRKNS_8DatatypeEii+0x9c>
    d479:	e8 52 d0 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d47e:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d482:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d486:	c9                   	leave  
    d487:	c3                   	ret    

000000000000d488 <_ZNK3MPI4Comm10Rsend_initEPKviRKNS_8DatatypeEii>:
    d488:	f3 0f 1e fa          	endbr64 
    d48c:	55                   	push   %rbp
    d48d:	48 89 e5             	mov    %rsp,%rbp
    d490:	53                   	push   %rbx
    d491:	48 83 ec 48          	sub    $0x48,%rsp
    d495:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d499:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d49d:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d4a1:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d4a4:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d4a8:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d4ac:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d4b3:	00 00 
    d4b5:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d4b9:	31 c0                	xor    %eax,%eax
    d4bb:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d4bf:	48 8b 58 08          	mov    0x8(%rax),%rbx
    d4c3:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d4c7:	48 89 c7             	mov    %rax,%rdi
    d4ca:	e8 d1 e3 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d4cf:	48 89 c7             	mov    %rax,%rdi
    d4d2:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d4d5:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d4d8:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d4dc:	48 83 ec 08          	sub    $0x8,%rsp
    d4e0:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d4e4:	51                   	push   %rcx
    d4e5:	49 89 d9             	mov    %rbx,%r9
    d4e8:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d4ec:	89 d1                	mov    %edx,%ecx
    d4ee:	48 89 fa             	mov    %rdi,%rdx
    d4f1:	48 89 c7             	mov    %rax,%rdi
    d4f4:	e8 07 cd ff ff       	call   a200 <MPI_Rsend_init@plt>
    d4f9:	48 83 c4 10          	add    $0x10,%rsp
    d4fd:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    d501:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d505:	48 89 d6             	mov    %rdx,%rsi
    d508:	48 89 c7             	mov    %rax,%rdi
    d50b:	e8 a6 e4 ff ff       	call   b9b6 <_ZN3MPI8PrequestC1ERKP14ompi_request_t>
    d510:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d514:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d51b:	00 00 
    d51d:	74 05                	je     d524 <_ZNK3MPI4Comm10Rsend_initEPKviRKNS_8DatatypeEii+0x9c>
    d51f:	e8 ac cf ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d524:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d528:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d52c:	c9                   	leave  
    d52d:	c3                   	ret    

000000000000d52e <_ZNK3MPI4Comm9Recv_initEPviRKNS_8DatatypeEii>:
    d52e:	f3 0f 1e fa          	endbr64 
    d532:	55                   	push   %rbp
    d533:	48 89 e5             	mov    %rsp,%rbp
    d536:	53                   	push   %rbx
    d537:	48 83 ec 48          	sub    $0x48,%rsp
    d53b:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d53f:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d543:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    d547:	89 4d c4             	mov    %ecx,-0x3c(%rbp)
    d54a:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    d54e:	44 89 4d c0          	mov    %r9d,-0x40(%rbp)
    d552:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d559:	00 00 
    d55b:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    d55f:	31 c0                	xor    %eax,%eax
    d561:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d565:	48 8b 58 08          	mov    0x8(%rax),%rbx
    d569:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    d56d:	48 89 c7             	mov    %rax,%rdi
    d570:	e8 2b e3 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d575:	48 89 c7             	mov    %rax,%rdi
    d578:	8b 55 c0             	mov    -0x40(%rbp),%edx
    d57b:	8b 75 c4             	mov    -0x3c(%rbp),%esi
    d57e:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    d582:	48 83 ec 08          	sub    $0x8,%rsp
    d586:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    d58a:	51                   	push   %rcx
    d58b:	49 89 d9             	mov    %rbx,%r9
    d58e:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    d592:	89 d1                	mov    %edx,%ecx
    d594:	48 89 fa             	mov    %rdi,%rdx
    d597:	48 89 c7             	mov    %rax,%rdi
    d59a:	e8 61 c8 ff ff       	call   9e00 <MPI_Recv_init@plt>
    d59f:	48 83 c4 10          	add    $0x10,%rsp
    d5a3:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    d5a7:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d5ab:	48 89 d6             	mov    %rdx,%rsi
    d5ae:	48 89 c7             	mov    %rax,%rdi
    d5b1:	e8 00 e4 ff ff       	call   b9b6 <_ZN3MPI8PrequestC1ERKP14ompi_request_t>
    d5b6:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d5ba:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d5c1:	00 00 
    d5c3:	74 05                	je     d5ca <_ZNK3MPI4Comm9Recv_initEPviRKNS_8DatatypeEii+0x9c>
    d5c5:	e8 06 cf ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d5ca:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d5ce:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d5d2:	c9                   	leave  
    d5d3:	c3                   	ret    

000000000000d5d4 <_ZNK3MPI4Comm8SendrecvEPKviRKNS_8DatatypeEiiPviS5_iiRNS_6StatusE>:
    d5d4:	f3 0f 1e fa          	endbr64 
    d5d8:	55                   	push   %rbp
    d5d9:	48 89 e5             	mov    %rsp,%rbp
    d5dc:	41 55                	push   %r13
    d5de:	41 54                	push   %r12
    d5e0:	53                   	push   %rbx
    d5e1:	48 83 ec 38          	sub    $0x38,%rsp
    d5e5:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    d5e9:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    d5ed:	89 55 cc             	mov    %edx,-0x34(%rbp)
    d5f0:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    d5f4:	44 89 45 c8          	mov    %r8d,-0x38(%rbp)
    d5f8:	44 89 4d bc          	mov    %r9d,-0x44(%rbp)
    d5fc:	48 8b 45 38          	mov    0x38(%rbp),%rax
    d600:	4c 8d 68 08          	lea    0x8(%rax),%r13
    d604:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    d608:	4c 8b 60 08          	mov    0x8(%rax),%r12
    d60c:	48 8b 45 20          	mov    0x20(%rbp),%rax
    d610:	48 89 c7             	mov    %rax,%rdi
    d613:	e8 88 e2 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d618:	48 89 c3             	mov    %rax,%rbx
    d61b:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    d61f:	48 89 c7             	mov    %rax,%rdi
    d622:	e8 79 e2 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d627:	48 89 c7             	mov    %rax,%rdi
    d62a:	44 8b 45 bc          	mov    -0x44(%rbp),%r8d
    d62e:	8b 55 c8             	mov    -0x38(%rbp),%edx
    d631:	8b 75 cc             	mov    -0x34(%rbp),%esi
    d634:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d638:	41 55                	push   %r13
    d63a:	41 54                	push   %r12
    d63c:	8b 4d 30             	mov    0x30(%rbp),%ecx
    d63f:	51                   	push   %rcx
    d640:	8b 4d 28             	mov    0x28(%rbp),%ecx
    d643:	51                   	push   %rcx
    d644:	53                   	push   %rbx
    d645:	8b 4d 18             	mov    0x18(%rbp),%ecx
    d648:	51                   	push   %rcx
    d649:	4c 8b 4d 10          	mov    0x10(%rbp),%r9
    d64d:	89 d1                	mov    %edx,%ecx
    d64f:	48 89 fa             	mov    %rdi,%rdx
    d652:	48 89 c7             	mov    %rax,%rdi
    d655:	e8 e6 cb ff ff       	call   a240 <MPI_Sendrecv@plt>
    d65a:	48 83 c4 30          	add    $0x30,%rsp
    d65e:	90                   	nop
    d65f:	48 8d 65 e8          	lea    -0x18(%rbp),%rsp
    d663:	5b                   	pop    %rbx
    d664:	41 5c                	pop    %r12
    d666:	41 5d                	pop    %r13
    d668:	5d                   	pop    %rbp
    d669:	c3                   	ret    

000000000000d66a <_ZNK3MPI4Comm8SendrecvEPKviRKNS_8DatatypeEiiPviS5_ii>:
    d66a:	f3 0f 1e fa          	endbr64 
    d66e:	55                   	push   %rbp
    d66f:	48 89 e5             	mov    %rsp,%rbp
    d672:	41 54                	push   %r12
    d674:	53                   	push   %rbx
    d675:	48 83 ec 30          	sub    $0x30,%rsp
    d679:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d67d:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    d681:	89 55 dc             	mov    %edx,-0x24(%rbp)
    d684:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    d688:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    d68c:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
    d690:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d694:	4c 8b 60 08          	mov    0x8(%rax),%r12
    d698:	48 8b 45 20          	mov    0x20(%rbp),%rax
    d69c:	48 89 c7             	mov    %rax,%rdi
    d69f:	e8 fc e1 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d6a4:	48 89 c3             	mov    %rax,%rbx
    d6a7:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d6ab:	48 89 c7             	mov    %rax,%rdi
    d6ae:	e8 ed e1 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d6b3:	48 89 c7             	mov    %rax,%rdi
    d6b6:	44 8b 45 cc          	mov    -0x34(%rbp),%r8d
    d6ba:	8b 55 d8             	mov    -0x28(%rbp),%edx
    d6bd:	8b 75 dc             	mov    -0x24(%rbp),%esi
    d6c0:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    d6c4:	6a 00                	push   $0x0
    d6c6:	41 54                	push   %r12
    d6c8:	8b 4d 30             	mov    0x30(%rbp),%ecx
    d6cb:	51                   	push   %rcx
    d6cc:	8b 4d 28             	mov    0x28(%rbp),%ecx
    d6cf:	51                   	push   %rcx
    d6d0:	53                   	push   %rbx
    d6d1:	8b 4d 18             	mov    0x18(%rbp),%ecx
    d6d4:	51                   	push   %rcx
    d6d5:	4c 8b 4d 10          	mov    0x10(%rbp),%r9
    d6d9:	89 d1                	mov    %edx,%ecx
    d6db:	48 89 fa             	mov    %rdi,%rdx
    d6de:	48 89 c7             	mov    %rax,%rdi
    d6e1:	e8 5a cb ff ff       	call   a240 <MPI_Sendrecv@plt>
    d6e6:	48 83 c4 30          	add    $0x30,%rsp
    d6ea:	90                   	nop
    d6eb:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    d6ef:	5b                   	pop    %rbx
    d6f0:	41 5c                	pop    %r12
    d6f2:	5d                   	pop    %rbp
    d6f3:	c3                   	ret    

000000000000d6f4 <_ZNK3MPI4Comm16Sendrecv_replaceEPviRKNS_8DatatypeEiiiiRNS_6StatusE>:
    d6f4:	f3 0f 1e fa          	endbr64 
    d6f8:	55                   	push   %rbp
    d6f9:	48 89 e5             	mov    %rsp,%rbp
    d6fc:	41 54                	push   %r12
    d6fe:	53                   	push   %rbx
    d6ff:	48 83 ec 30          	sub    $0x30,%rsp
    d703:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d707:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    d70b:	89 55 dc             	mov    %edx,-0x24(%rbp)
    d70e:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    d712:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    d716:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
    d71a:	48 8b 45 20          	mov    0x20(%rbp),%rax
    d71e:	4c 8d 60 08          	lea    0x8(%rax),%r12
    d722:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d726:	48 8b 58 08          	mov    0x8(%rax),%rbx
    d72a:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d72e:	48 89 c7             	mov    %rax,%rdi
    d731:	e8 6a e1 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d736:	48 89 c7             	mov    %rax,%rdi
    d739:	44 8b 45 cc          	mov    -0x34(%rbp),%r8d
    d73d:	8b 55 d8             	mov    -0x28(%rbp),%edx
    d740:	8b 75 dc             	mov    -0x24(%rbp),%esi
    d743:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    d747:	48 83 ec 08          	sub    $0x8,%rsp
    d74b:	41 54                	push   %r12
    d74d:	53                   	push   %rbx
    d74e:	8b 4d 18             	mov    0x18(%rbp),%ecx
    d751:	51                   	push   %rcx
    d752:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
    d756:	89 d1                	mov    %edx,%ecx
    d758:	48 89 fa             	mov    %rdi,%rdx
    d75b:	48 89 c7             	mov    %rax,%rdi
    d75e:	e8 ed ce ff ff       	call   a650 <MPI_Sendrecv_replace@plt>
    d763:	48 83 c4 20          	add    $0x20,%rsp
    d767:	90                   	nop
    d768:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    d76c:	5b                   	pop    %rbx
    d76d:	41 5c                	pop    %r12
    d76f:	5d                   	pop    %rbp
    d770:	c3                   	ret    
    d771:	90                   	nop

000000000000d772 <_ZNK3MPI4Comm16Sendrecv_replaceEPviRKNS_8DatatypeEiiii>:
    d772:	f3 0f 1e fa          	endbr64 
    d776:	55                   	push   %rbp
    d777:	48 89 e5             	mov    %rsp,%rbp
    d77a:	53                   	push   %rbx
    d77b:	48 83 ec 38          	sub    $0x38,%rsp
    d77f:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d783:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    d787:	89 55 dc             	mov    %edx,-0x24(%rbp)
    d78a:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    d78e:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    d792:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
    d796:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d79a:	48 8b 58 08          	mov    0x8(%rax),%rbx
    d79e:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d7a2:	48 89 c7             	mov    %rax,%rdi
    d7a5:	e8 f6 e0 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d7aa:	48 89 c7             	mov    %rax,%rdi
    d7ad:	44 8b 45 cc          	mov    -0x34(%rbp),%r8d
    d7b1:	8b 55 d8             	mov    -0x28(%rbp),%edx
    d7b4:	8b 75 dc             	mov    -0x24(%rbp),%esi
    d7b7:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    d7bb:	48 83 ec 08          	sub    $0x8,%rsp
    d7bf:	6a 00                	push   $0x0
    d7c1:	53                   	push   %rbx
    d7c2:	8b 4d 18             	mov    0x18(%rbp),%ecx
    d7c5:	51                   	push   %rcx
    d7c6:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
    d7ca:	89 d1                	mov    %edx,%ecx
    d7cc:	48 89 fa             	mov    %rdi,%rdx
    d7cf:	48 89 c7             	mov    %rax,%rdi
    d7d2:	e8 79 ce ff ff       	call   a650 <MPI_Sendrecv_replace@plt>
    d7d7:	48 83 c4 20          	add    $0x20,%rsp
    d7db:	90                   	nop
    d7dc:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d7e0:	c9                   	leave  
    d7e1:	c3                   	ret    

000000000000d7e2 <_ZNK3MPI4Comm9Get_groupEv>:
    d7e2:	f3 0f 1e fa          	endbr64 
    d7e6:	55                   	push   %rbp
    d7e7:	48 89 e5             	mov    %rsp,%rbp
    d7ea:	48 83 ec 20          	sub    $0x20,%rsp
    d7ee:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d7f2:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    d7f6:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d7fd:	00 00 
    d7ff:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    d803:	31 c0                	xor    %eax,%eax
    d805:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    d809:	48 8b 40 08          	mov    0x8(%rax),%rax
    d80d:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    d811:	48 89 d6             	mov    %rdx,%rsi
    d814:	48 89 c7             	mov    %rax,%rdi
    d817:	e8 d4 ca ff ff       	call   a2f0 <MPI_Comm_group@plt>
    d81c:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    d820:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d824:	48 89 d6             	mov    %rdx,%rsi
    d827:	48 89 c7             	mov    %rax,%rdi
    d82a:	e8 7f e2 ff ff       	call   baae <_ZN3MPI5GroupC1EP12ompi_group_t>
    d82f:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    d833:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    d83a:	00 00 
    d83c:	74 05                	je     d843 <_ZNK3MPI4Comm9Get_groupEv+0x61>
    d83e:	e8 8d cc ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d843:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d847:	c9                   	leave  
    d848:	c3                   	ret    
    d849:	90                   	nop

000000000000d84a <_ZNK3MPI4Comm8Get_sizeEv>:
    d84a:	f3 0f 1e fa          	endbr64 
    d84e:	55                   	push   %rbp
    d84f:	48 89 e5             	mov    %rsp,%rbp
    d852:	48 83 ec 20          	sub    $0x20,%rsp
    d856:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d85a:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d861:	00 00 
    d863:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    d867:	31 c0                	xor    %eax,%eax
    d869:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d86d:	48 8b 40 08          	mov    0x8(%rax),%rax
    d871:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    d875:	48 89 d6             	mov    %rdx,%rsi
    d878:	48 89 c7             	mov    %rax,%rdi
    d87b:	e8 60 cc ff ff       	call   a4e0 <MPI_Comm_size@plt>
    d880:	8b 45 f4             	mov    -0xc(%rbp),%eax
    d883:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    d887:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    d88e:	00 00 
    d890:	74 05                	je     d897 <_ZNK3MPI4Comm8Get_sizeEv+0x4d>
    d892:	e8 39 cc ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d897:	c9                   	leave  
    d898:	c3                   	ret    
    d899:	90                   	nop

000000000000d89a <_ZNK3MPI4Comm8Get_rankEv>:
    d89a:	f3 0f 1e fa          	endbr64 
    d89e:	55                   	push   %rbp
    d89f:	48 89 e5             	mov    %rsp,%rbp
    d8a2:	48 83 ec 20          	sub    $0x20,%rsp
    d8a6:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d8aa:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d8b1:	00 00 
    d8b3:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    d8b7:	31 c0                	xor    %eax,%eax
    d8b9:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d8bd:	48 8b 40 08          	mov    0x8(%rax),%rax
    d8c1:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    d8c5:	48 89 d6             	mov    %rdx,%rsi
    d8c8:	48 89 c7             	mov    %rax,%rdi
    d8cb:	e8 20 ce ff ff       	call   a6f0 <MPI_Comm_rank@plt>
    d8d0:	8b 45 f4             	mov    -0xc(%rbp),%eax
    d8d3:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    d8d7:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    d8de:	00 00 
    d8e0:	74 05                	je     d8e7 <_ZNK3MPI4Comm8Get_rankEv+0x4d>
    d8e2:	e8 e9 cb ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d8e7:	c9                   	leave  
    d8e8:	c3                   	ret    
    d8e9:	90                   	nop

000000000000d8ea <_ZN3MPI4Comm4FreeEv>:
    d8ea:	f3 0f 1e fa          	endbr64 
    d8ee:	55                   	push   %rbp
    d8ef:	48 89 e5             	mov    %rsp,%rbp
    d8f2:	48 83 ec 10          	sub    $0x10,%rsp
    d8f6:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    d8fa:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    d8fe:	48 83 c0 08          	add    $0x8,%rax
    d902:	48 89 c7             	mov    %rax,%rdi
    d905:	e8 36 cf ff ff       	call   a840 <MPI_Comm_free@plt>
    d90a:	90                   	nop
    d90b:	c9                   	leave  
    d90c:	c3                   	ret    
    d90d:	90                   	nop

000000000000d90e <_ZNK3MPI4Comm8Is_interEv>:
    d90e:	f3 0f 1e fa          	endbr64 
    d912:	55                   	push   %rbp
    d913:	48 89 e5             	mov    %rsp,%rbp
    d916:	48 83 ec 20          	sub    $0x20,%rsp
    d91a:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d91e:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    d925:	00 00 
    d927:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    d92b:	31 c0                	xor    %eax,%eax
    d92d:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d931:	48 8b 40 08          	mov    0x8(%rax),%rax
    d935:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    d939:	48 89 d6             	mov    %rdx,%rsi
    d93c:	48 89 c7             	mov    %rax,%rdi
    d93f:	e8 9c d1 ff ff       	call   aae0 <MPI_Comm_test_inter@plt>
    d944:	8b 45 f4             	mov    -0xc(%rbp),%eax
    d947:	85 c0                	test   %eax,%eax
    d949:	0f 95 c0             	setne  %al
    d94c:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    d950:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    d957:	00 00 
    d959:	74 05                	je     d960 <_ZNK3MPI4Comm8Is_interEv+0x52>
    d95b:	e8 70 cb ff ff       	call   a4d0 <__stack_chk_fail@plt>
    d960:	c9                   	leave  
    d961:	c3                   	ret    

000000000000d962 <_ZNK3MPI4Comm7BarrierEv>:
    d962:	f3 0f 1e fa          	endbr64 
    d966:	55                   	push   %rbp
    d967:	48 89 e5             	mov    %rsp,%rbp
    d96a:	48 83 ec 10          	sub    $0x10,%rsp
    d96e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    d972:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    d976:	48 8b 40 08          	mov    0x8(%rax),%rax
    d97a:	48 89 c7             	mov    %rax,%rdi
    d97d:	e8 3e c8 ff ff       	call   a1c0 <MPI_Barrier@plt>
    d982:	90                   	nop
    d983:	c9                   	leave  
    d984:	c3                   	ret    
    d985:	90                   	nop

000000000000d986 <_ZNK3MPI4Comm5BcastEPviRKNS_8DatatypeEi>:
    d986:	f3 0f 1e fa          	endbr64 
    d98a:	55                   	push   %rbp
    d98b:	48 89 e5             	mov    %rsp,%rbp
    d98e:	53                   	push   %rbx
    d98f:	48 83 ec 28          	sub    $0x28,%rsp
    d993:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d997:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    d99b:	89 55 dc             	mov    %edx,-0x24(%rbp)
    d99e:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    d9a2:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    d9a6:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    d9aa:	48 8b 58 08          	mov    0x8(%rax),%rbx
    d9ae:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    d9b2:	48 89 c7             	mov    %rax,%rdi
    d9b5:	e8 e6 de ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    d9ba:	48 89 c7             	mov    %rax,%rdi
    d9bd:	8b 55 d8             	mov    -0x28(%rbp),%edx
    d9c0:	8b 75 dc             	mov    -0x24(%rbp),%esi
    d9c3:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    d9c7:	49 89 d8             	mov    %rbx,%r8
    d9ca:	89 d1                	mov    %edx,%ecx
    d9cc:	48 89 fa             	mov    %rdi,%rdx
    d9cf:	48 89 c7             	mov    %rax,%rdi
    d9d2:	e8 79 c4 ff ff       	call   9e50 <MPI_Bcast@plt>
    d9d7:	90                   	nop
    d9d8:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    d9dc:	c9                   	leave  
    d9dd:	c3                   	ret    

000000000000d9de <_ZNK3MPI4Comm6GatherEPKviRKNS_8DatatypeEPviS5_i>:
    d9de:	f3 0f 1e fa          	endbr64 
    d9e2:	55                   	push   %rbp
    d9e3:	48 89 e5             	mov    %rsp,%rbp
    d9e6:	41 54                	push   %r12
    d9e8:	53                   	push   %rbx
    d9e9:	48 83 ec 30          	sub    $0x30,%rsp
    d9ed:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    d9f1:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    d9f5:	89 55 dc             	mov    %edx,-0x24(%rbp)
    d9f8:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    d9fc:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    da00:	44 89 4d d8          	mov    %r9d,-0x28(%rbp)
    da04:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    da08:	48 8b 58 08          	mov    0x8(%rax),%rbx
    da0c:	48 8b 7d 10          	mov    0x10(%rbp),%rdi
    da10:	e8 8b de ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    da15:	49 89 c4             	mov    %rax,%r12
    da18:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    da1c:	48 89 c7             	mov    %rax,%rdi
    da1f:	e8 7c de ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    da24:	48 89 c7             	mov    %rax,%rdi
    da27:	44 8b 45 d8          	mov    -0x28(%rbp),%r8d
    da2b:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    da2f:	8b 75 dc             	mov    -0x24(%rbp),%esi
    da32:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    da36:	53                   	push   %rbx
    da37:	8b 4d 18             	mov    0x18(%rbp),%ecx
    da3a:	51                   	push   %rcx
    da3b:	4d 89 e1             	mov    %r12,%r9
    da3e:	48 89 d1             	mov    %rdx,%rcx
    da41:	48 89 fa             	mov    %rdi,%rdx
    da44:	48 89 c7             	mov    %rax,%rdi
    da47:	e8 34 c8 ff ff       	call   a280 <MPI_Gather@plt>
    da4c:	48 83 c4 10          	add    $0x10,%rsp
    da50:	90                   	nop
    da51:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    da55:	5b                   	pop    %rbx
    da56:	41 5c                	pop    %r12
    da58:	5d                   	pop    %rbp
    da59:	c3                   	ret    

000000000000da5a <_ZNK3MPI4Comm7GathervEPKviRKNS_8DatatypeEPvPKiS8_S5_i>:
    da5a:	f3 0f 1e fa          	endbr64 
    da5e:	55                   	push   %rbp
    da5f:	48 89 e5             	mov    %rsp,%rbp
    da62:	41 54                	push   %r12
    da64:	53                   	push   %rbx
    da65:	48 83 ec 30          	sub    $0x30,%rsp
    da69:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    da6d:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    da71:	89 55 dc             	mov    %edx,-0x24(%rbp)
    da74:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    da78:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    da7c:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    da80:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    da84:	4c 8b 60 08          	mov    0x8(%rax),%r12
    da88:	48 8b 45 18          	mov    0x18(%rbp),%rax
    da8c:	48 89 c7             	mov    %rax,%rdi
    da8f:	e8 0c de ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    da94:	48 89 c3             	mov    %rax,%rbx
    da97:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    da9b:	48 89 c7             	mov    %rax,%rdi
    da9e:	e8 fd dd ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    daa3:	48 89 c7             	mov    %rax,%rdi
    daa6:	4c 8b 45 c0          	mov    -0x40(%rbp),%r8
    daaa:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    daae:	8b 75 dc             	mov    -0x24(%rbp),%esi
    dab1:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    dab5:	48 83 ec 08          	sub    $0x8,%rsp
    dab9:	41 54                	push   %r12
    dabb:	8b 4d 20             	mov    0x20(%rbp),%ecx
    dabe:	51                   	push   %rcx
    dabf:	53                   	push   %rbx
    dac0:	4c 8b 4d 10          	mov    0x10(%rbp),%r9
    dac4:	48 89 d1             	mov    %rdx,%rcx
    dac7:	48 89 fa             	mov    %rdi,%rdx
    daca:	48 89 c7             	mov    %rax,%rdi
    dacd:	e8 1e cf ff ff       	call   a9f0 <MPI_Gatherv@plt>
    dad2:	48 83 c4 20          	add    $0x20,%rsp
    dad6:	90                   	nop
    dad7:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    dadb:	5b                   	pop    %rbx
    dadc:	41 5c                	pop    %r12
    dade:	5d                   	pop    %rbp
    dadf:	c3                   	ret    

000000000000dae0 <_ZNK3MPI4Comm7ScatterEPKviRKNS_8DatatypeEPviS5_i>:
    dae0:	f3 0f 1e fa          	endbr64 
    dae4:	55                   	push   %rbp
    dae5:	48 89 e5             	mov    %rsp,%rbp
    dae8:	41 54                	push   %r12
    daea:	53                   	push   %rbx
    daeb:	48 83 ec 30          	sub    $0x30,%rsp
    daef:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    daf3:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    daf7:	89 55 dc             	mov    %edx,-0x24(%rbp)
    dafa:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    dafe:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    db02:	44 89 4d d8          	mov    %r9d,-0x28(%rbp)
    db06:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    db0a:	48 8b 58 08          	mov    0x8(%rax),%rbx
    db0e:	48 8b 7d 10          	mov    0x10(%rbp),%rdi
    db12:	e8 89 dd ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    db17:	49 89 c4             	mov    %rax,%r12
    db1a:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    db1e:	48 89 c7             	mov    %rax,%rdi
    db21:	e8 7a dd ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    db26:	48 89 c7             	mov    %rax,%rdi
    db29:	44 8b 45 d8          	mov    -0x28(%rbp),%r8d
    db2d:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    db31:	8b 75 dc             	mov    -0x24(%rbp),%esi
    db34:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    db38:	53                   	push   %rbx
    db39:	8b 4d 18             	mov    0x18(%rbp),%ecx
    db3c:	51                   	push   %rcx
    db3d:	4d 89 e1             	mov    %r12,%r9
    db40:	48 89 d1             	mov    %rdx,%rcx
    db43:	48 89 fa             	mov    %rdi,%rdx
    db46:	48 89 c7             	mov    %rax,%rdi
    db49:	e8 02 c6 ff ff       	call   a150 <MPI_Scatter@plt>
    db4e:	48 83 c4 10          	add    $0x10,%rsp
    db52:	90                   	nop
    db53:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    db57:	5b                   	pop    %rbx
    db58:	41 5c                	pop    %r12
    db5a:	5d                   	pop    %rbp
    db5b:	c3                   	ret    

000000000000db5c <_ZNK3MPI4Comm8ScattervEPKvPKiS4_RKNS_8DatatypeEPviS7_i>:
    db5c:	f3 0f 1e fa          	endbr64 
    db60:	55                   	push   %rbp
    db61:	48 89 e5             	mov    %rsp,%rbp
    db64:	41 54                	push   %r12
    db66:	53                   	push   %rbx
    db67:	48 83 ec 30          	sub    $0x30,%rsp
    db6b:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    db6f:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    db73:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    db77:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    db7b:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    db7f:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    db83:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    db87:	4c 8b 60 08          	mov    0x8(%rax),%r12
    db8b:	48 8b 45 18          	mov    0x18(%rbp),%rax
    db8f:	48 89 c7             	mov    %rax,%rdi
    db92:	e8 09 dd ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    db97:	48 89 c3             	mov    %rax,%rbx
    db9a:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    db9e:	48 89 c7             	mov    %rax,%rdi
    dba1:	e8 fa dc ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dba6:	48 89 c7             	mov    %rax,%rdi
    dba9:	4c 8b 45 c0          	mov    -0x40(%rbp),%r8
    dbad:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    dbb1:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    dbb5:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    dbb9:	48 83 ec 08          	sub    $0x8,%rsp
    dbbd:	41 54                	push   %r12
    dbbf:	8b 4d 20             	mov    0x20(%rbp),%ecx
    dbc2:	51                   	push   %rcx
    dbc3:	53                   	push   %rbx
    dbc4:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
    dbc8:	48 89 f9             	mov    %rdi,%rcx
    dbcb:	48 89 c7             	mov    %rax,%rdi
    dbce:	e8 9d c4 ff ff       	call   a070 <MPI_Scatterv@plt>
    dbd3:	48 83 c4 20          	add    $0x20,%rsp
    dbd7:	90                   	nop
    dbd8:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    dbdc:	5b                   	pop    %rbx
    dbdd:	41 5c                	pop    %r12
    dbdf:	5d                   	pop    %rbp
    dbe0:	c3                   	ret    
    dbe1:	90                   	nop

000000000000dbe2 <_ZNK3MPI4Comm9AllgatherEPKviRKNS_8DatatypeEPviS5_>:
    dbe2:	f3 0f 1e fa          	endbr64 
    dbe6:	55                   	push   %rbp
    dbe7:	48 89 e5             	mov    %rsp,%rbp
    dbea:	41 54                	push   %r12
    dbec:	53                   	push   %rbx
    dbed:	48 83 ec 30          	sub    $0x30,%rsp
    dbf1:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    dbf5:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    dbf9:	89 55 dc             	mov    %edx,-0x24(%rbp)
    dbfc:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    dc00:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    dc04:	44 89 4d d8          	mov    %r9d,-0x28(%rbp)
    dc08:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    dc0c:	48 8b 58 08          	mov    0x8(%rax),%rbx
    dc10:	48 8b 7d 10          	mov    0x10(%rbp),%rdi
    dc14:	e8 87 dc ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dc19:	49 89 c4             	mov    %rax,%r12
    dc1c:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    dc20:	48 89 c7             	mov    %rax,%rdi
    dc23:	e8 78 dc ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dc28:	48 89 c7             	mov    %rax,%rdi
    dc2b:	8b 4d d8             	mov    -0x28(%rbp),%ecx
    dc2e:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    dc32:	8b 75 dc             	mov    -0x24(%rbp),%esi
    dc35:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    dc39:	48 83 ec 08          	sub    $0x8,%rsp
    dc3d:	53                   	push   %rbx
    dc3e:	4d 89 e1             	mov    %r12,%r9
    dc41:	41 89 c8             	mov    %ecx,%r8d
    dc44:	48 89 d1             	mov    %rdx,%rcx
    dc47:	48 89 fa             	mov    %rdi,%rdx
    dc4a:	48 89 c7             	mov    %rax,%rdi
    dc4d:	e8 9e c1 ff ff       	call   9df0 <MPI_Allgather@plt>
    dc52:	48 83 c4 10          	add    $0x10,%rsp
    dc56:	90                   	nop
    dc57:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    dc5b:	5b                   	pop    %rbx
    dc5c:	41 5c                	pop    %r12
    dc5e:	5d                   	pop    %rbp
    dc5f:	c3                   	ret    

000000000000dc60 <_ZNK3MPI4Comm10AllgathervEPKviRKNS_8DatatypeEPvPKiS8_S5_>:
    dc60:	f3 0f 1e fa          	endbr64 
    dc64:	55                   	push   %rbp
    dc65:	48 89 e5             	mov    %rsp,%rbp
    dc68:	41 54                	push   %r12
    dc6a:	53                   	push   %rbx
    dc6b:	48 83 ec 30          	sub    $0x30,%rsp
    dc6f:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    dc73:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    dc77:	89 55 dc             	mov    %edx,-0x24(%rbp)
    dc7a:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    dc7e:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    dc82:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    dc86:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    dc8a:	4c 8b 60 08          	mov    0x8(%rax),%r12
    dc8e:	48 8b 45 18          	mov    0x18(%rbp),%rax
    dc92:	48 89 c7             	mov    %rax,%rdi
    dc95:	e8 06 dc ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dc9a:	48 89 c3             	mov    %rax,%rbx
    dc9d:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    dca1:	48 89 c7             	mov    %rax,%rdi
    dca4:	e8 f7 db ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dca9:	48 89 c7             	mov    %rax,%rdi
    dcac:	48 8b 4d c0          	mov    -0x40(%rbp),%rcx
    dcb0:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    dcb4:	8b 75 dc             	mov    -0x24(%rbp),%esi
    dcb7:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    dcbb:	41 54                	push   %r12
    dcbd:	53                   	push   %rbx
    dcbe:	4c 8b 4d 10          	mov    0x10(%rbp),%r9
    dcc2:	49 89 c8             	mov    %rcx,%r8
    dcc5:	48 89 d1             	mov    %rdx,%rcx
    dcc8:	48 89 fa             	mov    %rdi,%rdx
    dccb:	48 89 c7             	mov    %rax,%rdi
    dcce:	e8 8d c8 ff ff       	call   a560 <MPI_Allgatherv@plt>
    dcd3:	48 83 c4 10          	add    $0x10,%rsp
    dcd7:	90                   	nop
    dcd8:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    dcdc:	5b                   	pop    %rbx
    dcdd:	41 5c                	pop    %r12
    dcdf:	5d                   	pop    %rbp
    dce0:	c3                   	ret    
    dce1:	90                   	nop

000000000000dce2 <_ZNK3MPI4Comm8AlltoallEPKviRKNS_8DatatypeEPviS5_>:
    dce2:	f3 0f 1e fa          	endbr64 
    dce6:	55                   	push   %rbp
    dce7:	48 89 e5             	mov    %rsp,%rbp
    dcea:	41 54                	push   %r12
    dcec:	53                   	push   %rbx
    dced:	48 83 ec 30          	sub    $0x30,%rsp
    dcf1:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    dcf5:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    dcf9:	89 55 dc             	mov    %edx,-0x24(%rbp)
    dcfc:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    dd00:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    dd04:	44 89 4d d8          	mov    %r9d,-0x28(%rbp)
    dd08:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    dd0c:	48 8b 58 08          	mov    0x8(%rax),%rbx
    dd10:	48 8b 7d 10          	mov    0x10(%rbp),%rdi
    dd14:	e8 87 db ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dd19:	49 89 c4             	mov    %rax,%r12
    dd1c:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    dd20:	48 89 c7             	mov    %rax,%rdi
    dd23:	e8 78 db ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dd28:	48 89 c7             	mov    %rax,%rdi
    dd2b:	8b 4d d8             	mov    -0x28(%rbp),%ecx
    dd2e:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    dd32:	8b 75 dc             	mov    -0x24(%rbp),%esi
    dd35:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    dd39:	48 83 ec 08          	sub    $0x8,%rsp
    dd3d:	53                   	push   %rbx
    dd3e:	4d 89 e1             	mov    %r12,%r9
    dd41:	41 89 c8             	mov    %ecx,%r8d
    dd44:	48 89 d1             	mov    %rdx,%rcx
    dd47:	48 89 fa             	mov    %rdi,%rdx
    dd4a:	48 89 c7             	mov    %rax,%rdi
    dd4d:	e8 ee c2 ff ff       	call   a040 <MPI_Alltoall@plt>
    dd52:	48 83 c4 10          	add    $0x10,%rsp
    dd56:	90                   	nop
    dd57:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    dd5b:	5b                   	pop    %rbx
    dd5c:	41 5c                	pop    %r12
    dd5e:	5d                   	pop    %rbp
    dd5f:	c3                   	ret    

000000000000dd60 <_ZNK3MPI4Comm9AlltoallvEPKvPKiS4_RKNS_8DatatypeEPvS4_S4_S7_>:
    dd60:	f3 0f 1e fa          	endbr64 
    dd64:	55                   	push   %rbp
    dd65:	48 89 e5             	mov    %rsp,%rbp
    dd68:	41 54                	push   %r12
    dd6a:	53                   	push   %rbx
    dd6b:	48 83 ec 30          	sub    $0x30,%rsp
    dd6f:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    dd73:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    dd77:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    dd7b:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    dd7f:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    dd83:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    dd87:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    dd8b:	4c 8b 60 08          	mov    0x8(%rax),%r12
    dd8f:	48 8b 45 20          	mov    0x20(%rbp),%rax
    dd93:	48 89 c7             	mov    %rax,%rdi
    dd96:	e8 05 db ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    dd9b:	48 89 c3             	mov    %rax,%rbx
    dd9e:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    dda2:	48 89 c7             	mov    %rax,%rdi
    dda5:	e8 f6 da ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    ddaa:	48 89 c1             	mov    %rax,%rcx
    ddad:	48 8b 7d c0          	mov    -0x40(%rbp),%rdi
    ddb1:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    ddb5:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    ddb9:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    ddbd:	48 83 ec 08          	sub    $0x8,%rsp
    ddc1:	41 54                	push   %r12
    ddc3:	53                   	push   %rbx
    ddc4:	ff 75 18             	push   0x18(%rbp)
    ddc7:	4c 8b 4d 10          	mov    0x10(%rbp),%r9
    ddcb:	49 89 f8             	mov    %rdi,%r8
    ddce:	48 89 c7             	mov    %rax,%rdi
    ddd1:	e8 1a ca ff ff       	call   a7f0 <MPI_Alltoallv@plt>
    ddd6:	48 83 c4 20          	add    $0x20,%rsp
    ddda:	90                   	nop
    dddb:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    dddf:	5b                   	pop    %rbx
    dde0:	41 5c                	pop    %r12
    dde2:	5d                   	pop    %rbp
    dde3:	c3                   	ret    

000000000000dde4 <_ZNK3MPI4Comm9AlltoallwEPKvPKiS4_PKNS_8DatatypeEPvS4_S4_S7_>:
    dde4:	f3 0f 1e fa          	endbr64 
    dde8:	55                   	push   %rbp
    dde9:	48 89 e5             	mov    %rsp,%rbp
    ddec:	53                   	push   %rbx
    dded:	48 83 ec 48          	sub    $0x48,%rsp
    ddf1:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    ddf5:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    ddf9:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    ddfd:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    de01:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    de05:	4c 89 4d b0          	mov    %r9,-0x50(%rbp)
    de09:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    de0d:	48 8b 00             	mov    (%rax),%rax
    de10:	48 05 d8 00 00 00    	add    $0xd8,%rax
    de16:	48 8b 10             	mov    (%rax),%rdx
    de19:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    de1d:	48 89 c7             	mov    %rax,%rdi
    de20:	ff d2                	call   *%rdx
    de22:	89 45 e4             	mov    %eax,-0x1c(%rbp)
    de25:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    de28:	01 c0                	add    %eax,%eax
    de2a:	48 98                	cltq   
    de2c:	48 ba ff ff ff ff ff 	movabs $0xfffffffffffffff,%rdx
    de33:	ff ff 0f 
    de36:	48 39 d0             	cmp    %rdx,%rax
    de39:	77 20                	ja     de5b <_ZNK3MPI4Comm9AlltoallwEPKvPKiS4_PKNS_8DatatypeEPvS4_S4_S7_+0x77>
    de3b:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    de3e:	01 c0                	add    %eax,%eax
    de40:	48 98                	cltq   
    de42:	48 c1 e0 03          	shl    $0x3,%rax
    de46:	48 89 c7             	mov    %rax,%rdi
    de49:	e8 52 bf ff ff       	call   9da0 <_Znam@plt>
    de4e:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    de52:	c7 45 e0 00 00 00 00 	movl   $0x0,-0x20(%rbp)
    de59:	eb 76                	jmp    ded1 <_ZNK3MPI4Comm9AlltoallwEPKvPKiS4_PKNS_8DatatypeEPvS4_S4_S7_+0xed>
    de5b:	e8 80 c8 ff ff       	call   a6e0 <__cxa_throw_bad_array_new_length@plt>
    de60:	8b 45 e0             	mov    -0x20(%rbp),%eax
    de63:	48 98                	cltq   
    de65:	48 c1 e0 04          	shl    $0x4,%rax
    de69:	48 89 c2             	mov    %rax,%rdx
    de6c:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    de70:	48 01 d0             	add    %rdx,%rax
    de73:	8b 55 e0             	mov    -0x20(%rbp),%edx
    de76:	48 63 d2             	movslq %edx,%rdx
    de79:	48 8d 0c d5 00 00 00 	lea    0x0(,%rdx,8),%rcx
    de80:	00 
    de81:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    de85:	48 8d 1c 11          	lea    (%rcx,%rdx,1),%rbx
    de89:	48 89 c7             	mov    %rax,%rdi
    de8c:	e8 0f da ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    de91:	48 89 03             	mov    %rax,(%rbx)
    de94:	8b 45 e0             	mov    -0x20(%rbp),%eax
    de97:	48 98                	cltq   
    de99:	48 c1 e0 04          	shl    $0x4,%rax
    de9d:	48 89 c2             	mov    %rax,%rdx
    dea0:	48 8b 45 20          	mov    0x20(%rbp),%rax
    dea4:	48 01 d0             	add    %rdx,%rax
    dea7:	8b 4d e0             	mov    -0x20(%rbp),%ecx
    deaa:	8b 55 e4             	mov    -0x1c(%rbp),%edx
    dead:	01 ca                	add    %ecx,%edx
    deaf:	48 63 d2             	movslq %edx,%rdx
    deb2:	48 8d 0c d5 00 00 00 	lea    0x0(,%rdx,8),%rcx
    deb9:	00 
    deba:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    debe:	48 8d 1c 11          	lea    (%rcx,%rdx,1),%rbx
    dec2:	48 89 c7             	mov    %rax,%rdi
    dec5:	e8 d6 d9 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    deca:	48 89 03             	mov    %rax,(%rbx)
    decd:	83 45 e0 01          	addl   $0x1,-0x20(%rbp)
    ded1:	8b 45 e0             	mov    -0x20(%rbp),%eax
    ded4:	3b 45 e4             	cmp    -0x1c(%rbp),%eax
    ded7:	7c 87                	jl     de60 <_ZNK3MPI4Comm9AlltoallwEPKvPKiS4_PKNS_8DatatypeEPvS4_S4_S7_+0x7c>
    ded9:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    dedd:	48 8b 78 08          	mov    0x8(%rax),%rdi
    dee1:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    dee4:	48 98                	cltq   
    dee6:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
    deed:	00 
    deee:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    def2:	4c 8d 04 02          	lea    (%rdx,%rax,1),%r8
    def6:	4c 8b 55 b0          	mov    -0x50(%rbp),%r10
    defa:	48 8b 4d e8          	mov    -0x18(%rbp),%rcx
    defe:	48 8b 55 c0          	mov    -0x40(%rbp),%rdx
    df02:	48 8b 75 c8          	mov    -0x38(%rbp),%rsi
    df06:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    df0a:	48 83 ec 08          	sub    $0x8,%rsp
    df0e:	57                   	push   %rdi
    df0f:	41 50                	push   %r8
    df11:	ff 75 18             	push   0x18(%rbp)
    df14:	4c 8b 4d 10          	mov    0x10(%rbp),%r9
    df18:	4d 89 d0             	mov    %r10,%r8
    df1b:	48 89 c7             	mov    %rax,%rdi
    df1e:	e8 5d c7 ff ff       	call   a680 <MPI_Alltoallw@plt>
    df23:	48 83 c4 20          	add    $0x20,%rsp
    df27:	48 83 7d e8 00       	cmpq   $0x0,-0x18(%rbp)
    df2c:	74 0c                	je     df3a <_ZNK3MPI4Comm9AlltoallwEPKvPKiS4_PKNS_8DatatypeEPvS4_S4_S7_+0x156>
    df2e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    df32:	48 89 c7             	mov    %rax,%rdi
    df35:	e8 86 c6 ff ff       	call   a5c0 <_ZdaPv@plt>
    df3a:	90                   	nop
    df3b:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    df3f:	c9                   	leave  
    df40:	c3                   	ret    
    df41:	90                   	nop

000000000000df42 <_ZNK3MPI4Comm6ReduceEPKvPviRKNS_8DatatypeERKNS_2OpEi>:
    df42:	f3 0f 1e fa          	endbr64 
    df46:	55                   	push   %rbp
    df47:	48 89 e5             	mov    %rsp,%rbp
    df4a:	41 54                	push   %r12
    df4c:	53                   	push   %rbx
    df4d:	48 83 ec 30          	sub    $0x30,%rsp
    df51:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    df55:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    df59:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    df5d:	89 4d d4             	mov    %ecx,-0x2c(%rbp)
    df60:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    df64:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    df68:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    df6c:	48 8b 58 08          	mov    0x8(%rax),%rbx
    df70:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    df74:	48 89 c7             	mov    %rax,%rdi
    df77:	e8 9a 1b 00 00       	call   fb16 <_ZNK3MPI2OpcvP9ompi_op_tEv>
    df7c:	49 89 c4             	mov    %rax,%r12
    df7f:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    df83:	48 89 c7             	mov    %rax,%rdi
    df86:	e8 15 d9 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    df8b:	48 89 c1             	mov    %rax,%rcx
    df8e:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    df91:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    df95:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    df99:	48 83 ec 08          	sub    $0x8,%rsp
    df9d:	53                   	push   %rbx
    df9e:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
    dfa2:	4d 89 e0             	mov    %r12,%r8
    dfa5:	48 89 c7             	mov    %rax,%rdi
    dfa8:	e8 53 c6 ff ff       	call   a600 <MPI_Reduce@plt>
    dfad:	48 83 c4 10          	add    $0x10,%rsp
    dfb1:	90                   	nop
    dfb2:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
    dfb6:	5b                   	pop    %rbx
    dfb7:	41 5c                	pop    %r12
    dfb9:	5d                   	pop    %rbp
    dfba:	c3                   	ret    
    dfbb:	90                   	nop

000000000000dfbc <_ZNK3MPI4Comm9AllreduceEPKvPviRKNS_8DatatypeERKNS_2OpE>:
    dfbc:	f3 0f 1e fa          	endbr64 
    dfc0:	55                   	push   %rbp
    dfc1:	48 89 e5             	mov    %rsp,%rbp
    dfc4:	41 54                	push   %r12
    dfc6:	53                   	push   %rbx
    dfc7:	48 83 ec 30          	sub    $0x30,%rsp
    dfcb:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    dfcf:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    dfd3:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    dfd7:	89 4d d4             	mov    %ecx,-0x2c(%rbp)
    dfda:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    dfde:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    dfe2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    dfe6:	4c 8b 60 08          	mov    0x8(%rax),%r12
    dfea:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    dfee:	48 89 c7             	mov    %rax,%rdi
    dff1:	e8 20 1b 00 00       	call   fb16 <_ZNK3MPI2OpcvP9ompi_op_tEv>
    dff6:	48 89 c3             	mov    %rax,%rbx
    dff9:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    dffd:	48 89 c7             	mov    %rax,%rdi
    e000:	e8 9b d8 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    e005:	48 89 c1             	mov    %rax,%rcx
    e008:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    e00b:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    e00f:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e013:	4d 89 e1             	mov    %r12,%r9
    e016:	49 89 d8             	mov    %rbx,%r8
    e019:	48 89 c7             	mov    %rax,%rdi
    e01c:	e8 2f c5 ff ff       	call   a550 <MPI_Allreduce@plt>
    e021:	90                   	nop
    e022:	48 83 c4 30          	add    $0x30,%rsp
    e026:	5b                   	pop    %rbx
    e027:	41 5c                	pop    %r12
    e029:	5d                   	pop    %rbp
    e02a:	c3                   	ret    
    e02b:	90                   	nop

000000000000e02c <_ZNK3MPI4Comm14Reduce_scatterEPKvPvPiRKNS_8DatatypeERKNS_2OpE>:
    e02c:	f3 0f 1e fa          	endbr64 
    e030:	55                   	push   %rbp
    e031:	48 89 e5             	mov    %rsp,%rbp
    e034:	41 54                	push   %r12
    e036:	53                   	push   %rbx
    e037:	48 83 ec 30          	sub    $0x30,%rsp
    e03b:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e03f:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e043:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    e047:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    e04b:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    e04f:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    e053:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e057:	4c 8b 60 08          	mov    0x8(%rax),%r12
    e05b:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    e05f:	48 89 c7             	mov    %rax,%rdi
    e062:	e8 af 1a 00 00       	call   fb16 <_ZNK3MPI2OpcvP9ompi_op_tEv>
    e067:	48 89 c3             	mov    %rax,%rbx
    e06a:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e06e:	48 89 c7             	mov    %rax,%rdi
    e071:	e8 2a d8 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    e076:	48 89 c1             	mov    %rax,%rcx
    e079:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    e07d:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    e081:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e085:	4d 89 e1             	mov    %r12,%r9
    e088:	49 89 d8             	mov    %rbx,%r8
    e08b:	48 89 c7             	mov    %rax,%rdi
    e08e:	e8 8d c4 ff ff       	call   a520 <MPI_Reduce_scatter@plt>
    e093:	90                   	nop
    e094:	48 83 c4 30          	add    $0x30,%rsp
    e098:	5b                   	pop    %rbx
    e099:	41 5c                	pop    %r12
    e09b:	5d                   	pop    %rbp
    e09c:	c3                   	ret    
    e09d:	90                   	nop

000000000000e09e <_ZN3MPI4Comm10DisconnectEv>:
    e09e:	f3 0f 1e fa          	endbr64 
    e0a2:	55                   	push   %rbp
    e0a3:	48 89 e5             	mov    %rsp,%rbp
    e0a6:	48 83 ec 10          	sub    $0x10,%rsp
    e0aa:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    e0ae:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e0b2:	48 83 c0 08          	add    $0x8,%rax
    e0b6:	48 89 c7             	mov    %rax,%rdi
    e0b9:	e8 02 be ff ff       	call   9ec0 <MPI_Comm_disconnect@plt>
    e0be:	90                   	nop
    e0bf:	c9                   	leave  
    e0c0:	c3                   	ret    
    e0c1:	90                   	nop

000000000000e0c2 <_ZNK3MPI4Comm8Get_nameEPcRi>:
    e0c2:	f3 0f 1e fa          	endbr64 
    e0c6:	55                   	push   %rbp
    e0c7:	48 89 e5             	mov    %rsp,%rbp
    e0ca:	48 83 ec 20          	sub    $0x20,%rsp
    e0ce:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    e0d2:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    e0d6:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    e0da:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e0de:	48 8b 40 08          	mov    0x8(%rax),%rax
    e0e2:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    e0e6:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    e0ea:	48 89 ce             	mov    %rcx,%rsi
    e0ed:	48 89 c7             	mov    %rax,%rdi
    e0f0:	e8 2b bd ff ff       	call   9e20 <MPI_Comm_get_name@plt>
    e0f5:	90                   	nop
    e0f6:	c9                   	leave  
    e0f7:	c3                   	ret    

000000000000e0f8 <_ZN3MPI4Comm8Set_nameEPKc>:
    e0f8:	f3 0f 1e fa          	endbr64 
    e0fc:	55                   	push   %rbp
    e0fd:	48 89 e5             	mov    %rsp,%rbp
    e100:	48 83 ec 10          	sub    $0x10,%rsp
    e104:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    e108:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    e10c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e110:	48 8b 40 08          	mov    0x8(%rax),%rax
    e114:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    e118:	48 89 d6             	mov    %rdx,%rsi
    e11b:	48 89 c7             	mov    %rax,%rdi
    e11e:	e8 fd c7 ff ff       	call   a920 <MPI_Comm_set_name@plt>
    e123:	90                   	nop
    e124:	c9                   	leave  
    e125:	c3                   	ret    

000000000000e126 <_ZNK3MPI4Comm12Get_topologyEv>:
    e126:	f3 0f 1e fa          	endbr64 
    e12a:	55                   	push   %rbp
    e12b:	48 89 e5             	mov    %rsp,%rbp
    e12e:	48 83 ec 20          	sub    $0x20,%rsp
    e132:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e136:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e13d:	00 00 
    e13f:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e143:	31 c0                	xor    %eax,%eax
    e145:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e149:	48 8b 40 08          	mov    0x8(%rax),%rax
    e14d:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    e151:	48 89 d6             	mov    %rdx,%rsi
    e154:	48 89 c7             	mov    %rax,%rdi
    e157:	e8 d4 c5 ff ff       	call   a730 <MPI_Topo_test@plt>
    e15c:	8b 45 f4             	mov    -0xc(%rbp),%eax
    e15f:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    e163:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    e16a:	00 00 
    e16c:	74 05                	je     e173 <_ZNK3MPI4Comm12Get_topologyEv+0x4d>
    e16e:	e8 5d c3 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    e173:	c9                   	leave  
    e174:	c3                   	ret    
    e175:	90                   	nop

000000000000e176 <_ZN3MPI4Comm5AbortEi>:
    e176:	f3 0f 1e fa          	endbr64 
    e17a:	55                   	push   %rbp
    e17b:	48 89 e5             	mov    %rsp,%rbp
    e17e:	48 83 ec 10          	sub    $0x10,%rsp
    e182:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    e186:	89 75 f4             	mov    %esi,-0xc(%rbp)
    e189:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e18d:	48 8b 40 08          	mov    0x8(%rax),%rax
    e191:	8b 55 f4             	mov    -0xc(%rbp),%edx
    e194:	89 d6                	mov    %edx,%esi
    e196:	48 89 c7             	mov    %rax,%rdi
    e199:	e8 72 c2 ff ff       	call   a410 <MPI_Abort@plt>
    e19e:	90                   	nop
    e19f:	c9                   	leave  
    e1a0:	c3                   	ret    
    e1a1:	90                   	nop

000000000000e1a2 <_ZNK3MPI4Comm14Get_errhandlerEv>:
    e1a2:	f3 0f 1e fa          	endbr64 
    e1a6:	55                   	push   %rbp
    e1a7:	48 89 e5             	mov    %rsp,%rbp
    e1aa:	48 83 ec 20          	sub    $0x20,%rsp
    e1ae:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e1b2:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e1b6:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e1bd:	00 00 
    e1bf:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e1c3:	31 c0                	xor    %eax,%eax
    e1c5:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e1c9:	48 8b 40 08          	mov    0x8(%rax),%rax
    e1cd:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    e1d1:	48 89 d6             	mov    %rdx,%rsi
    e1d4:	48 89 c7             	mov    %rax,%rdi
    e1d7:	e8 14 c4 ff ff       	call   a5f0 <MPI_Comm_get_errhandler@plt>
    e1dc:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    e1e0:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e1e4:	48 89 d6             	mov    %rdx,%rsi
    e1e7:	48 89 c7             	mov    %rax,%rdi
    e1ea:	e8 b7 da ff ff       	call   bca6 <_ZN3MPI10ErrhandlerC1EP17ompi_errhandler_t>
    e1ef:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e1f3:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e1fa:	00 00 
    e1fc:	74 05                	je     e203 <_ZNK3MPI4Comm14Get_errhandlerEv+0x61>
    e1fe:	e8 cd c2 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    e203:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e207:	c9                   	leave  
    e208:	c3                   	ret    
    e209:	90                   	nop

000000000000e20a <_ZN3MPI4Comm14Set_errhandlerERKNS_10ErrhandlerE>:
    e20a:	f3 0f 1e fa          	endbr64 
    e20e:	55                   	push   %rbp
    e20f:	48 89 e5             	mov    %rsp,%rbp
    e212:	48 83 ec 10          	sub    $0x10,%rsp
    e216:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    e21a:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    e21e:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    e222:	48 89 c7             	mov    %rax,%rdi
    e225:	e8 aa da ff ff       	call   bcd4 <_ZNK3MPI10ErrhandlercvP17ompi_errhandler_tEv>
    e22a:	48 89 c2             	mov    %rax,%rdx
    e22d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e231:	48 8b 40 08          	mov    0x8(%rax),%rax
    e235:	48 89 d6             	mov    %rdx,%rsi
    e238:	48 89 c7             	mov    %rax,%rdi
    e23b:	e8 90 c5 ff ff       	call   a7d0 <MPI_Comm_set_errhandler@plt>
    e240:	90                   	nop
    e241:	c9                   	leave  
    e242:	c3                   	ret    
    e243:	90                   	nop

000000000000e244 <_ZNK3MPI4Comm8Set_attrEiPKv>:
    e244:	f3 0f 1e fa          	endbr64 
    e248:	55                   	push   %rbp
    e249:	48 89 e5             	mov    %rsp,%rbp
    e24c:	48 83 ec 20          	sub    $0x20,%rsp
    e250:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    e254:	89 75 f4             	mov    %esi,-0xc(%rbp)
    e257:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    e25b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e25f:	48 8b 40 08          	mov    0x8(%rax),%rax
    e263:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    e267:	8b 4d f4             	mov    -0xc(%rbp),%ecx
    e26a:	89 ce                	mov    %ecx,%esi
    e26c:	48 89 c7             	mov    %rax,%rdi
    e26f:	e8 ac bd ff ff       	call   a020 <MPI_Comm_set_attr@plt>
    e274:	90                   	nop
    e275:	c9                   	leave  
    e276:	c3                   	ret    
    e277:	90                   	nop

000000000000e278 <_ZNK3MPI4Comm8Get_attrEiPv>:
    e278:	f3 0f 1e fa          	endbr64 
    e27c:	55                   	push   %rbp
    e27d:	48 89 e5             	mov    %rsp,%rbp
    e280:	48 83 ec 30          	sub    $0x30,%rsp
    e284:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e288:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    e28b:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    e28f:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e296:	00 00 
    e298:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e29c:	31 c0                	xor    %eax,%eax
    e29e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e2a2:	48 8b 40 08          	mov    0x8(%rax),%rax
    e2a6:	48 8d 4d f4          	lea    -0xc(%rbp),%rcx
    e2aa:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
    e2ae:	8b 75 e4             	mov    -0x1c(%rbp),%esi
    e2b1:	48 89 c7             	mov    %rax,%rdi
    e2b4:	e8 17 c4 ff ff       	call   a6d0 <MPI_Comm_get_attr@plt>
    e2b9:	8b 45 f4             	mov    -0xc(%rbp),%eax
    e2bc:	85 c0                	test   %eax,%eax
    e2be:	0f 95 c0             	setne  %al
    e2c1:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    e2c5:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    e2cc:	00 00 
    e2ce:	74 05                	je     e2d5 <_ZNK3MPI4Comm8Get_attrEiPv+0x5d>
    e2d0:	e8 fb c1 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    e2d5:	c9                   	leave  
    e2d6:	c3                   	ret    
    e2d7:	90                   	nop

000000000000e2d8 <_ZN3MPI4Comm11Delete_attrEi>:
    e2d8:	f3 0f 1e fa          	endbr64 
    e2dc:	55                   	push   %rbp
    e2dd:	48 89 e5             	mov    %rsp,%rbp
    e2e0:	48 83 ec 10          	sub    $0x10,%rsp
    e2e4:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    e2e8:	89 75 f4             	mov    %esi,-0xc(%rbp)
    e2eb:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e2ef:	48 8b 40 08          	mov    0x8(%rax),%rax
    e2f3:	8b 55 f4             	mov    -0xc(%rbp),%edx
    e2f6:	89 d6                	mov    %edx,%esi
    e2f8:	48 89 c7             	mov    %rax,%rdi
    e2fb:	e8 80 c1 ff ff       	call   a480 <MPI_Comm_delete_attr@plt>
    e300:	90                   	nop
    e301:	c9                   	leave  
    e302:	c3                   	ret    
    e303:	90                   	nop

000000000000e304 <_ZN3MPI9IntracommC1EP19ompi_communicator_t>:
    e304:	f3 0f 1e fa          	endbr64 
    e308:	55                   	push   %rbp
    e309:	48 89 e5             	mov    %rsp,%rbp
    e30c:	53                   	push   %rbx
    e30d:	48 83 ec 28          	sub    $0x28,%rsp
    e311:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    e315:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    e319:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e320:	00 00 
    e322:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    e326:	31 c0                	xor    %eax,%eax
    e328:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e32c:	48 89 c7             	mov    %rax,%rdi
    e32f:	e8 5c bd ff ff       	call   a090 <_ZN3MPI4CommC2Ev@plt>
    e334:	48 8d 15 2d db 00 00 	lea    0xdb2d(%rip),%rdx        # 1be68 <_ZTVN3MPI9IntracommE+0x10>
    e33b:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e33f:	48 89 10             	mov    %rdx,(%rax)
    e342:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
    e349:	e8 a0 e5 ff ff       	call   c8ee <_ZN3MPI14Is_initializedEv>
    e34e:	84 c0                	test   %al,%al
    e350:	74 14                	je     e366 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0x62>
    e352:	48 8d 05 c7 fa 00 00 	lea    0xfac7(%rip),%rax        # 1de20 <ompi_mpi_comm_null@@Base>
    e359:	48 39 45 d0          	cmp    %rax,-0x30(%rbp)
    e35d:	74 07                	je     e366 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0x62>
    e35f:	b8 01 00 00 00       	mov    $0x1,%eax
    e364:	eb 05                	jmp    e36b <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0x67>
    e366:	b8 00 00 00 00       	mov    $0x0,%eax
    e36b:	84 c0                	test   %al,%al
    e36d:	74 39                	je     e3a8 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0xa4>
    e36f:	48 8d 55 e4          	lea    -0x1c(%rbp),%rdx
    e373:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    e377:	48 89 d6             	mov    %rdx,%rsi
    e37a:	48 89 c7             	mov    %rax,%rdi
    e37d:	e8 5e c7 ff ff       	call   aae0 <MPI_Comm_test_inter@plt>
    e382:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    e385:	85 c0                	test   %eax,%eax
    e387:	74 11                	je     e39a <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0x96>
    e389:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e38d:	48 8d 15 8c fa 00 00 	lea    0xfa8c(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    e394:	48 89 50 08          	mov    %rdx,0x8(%rax)
    e398:	eb 3a                	jmp    e3d4 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0xd0>
    e39a:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e39e:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    e3a2:	48 89 50 08          	mov    %rdx,0x8(%rax)
    e3a6:	eb 2c                	jmp    e3d4 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0xd0>
    e3a8:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e3ac:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    e3b0:	48 89 50 08          	mov    %rdx,0x8(%rax)
    e3b4:	eb 1e                	jmp    e3d4 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0xd0>
    e3b6:	f3 0f 1e fa          	endbr64 
    e3ba:	48 89 c3             	mov    %rax,%rbx
    e3bd:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e3c1:	48 89 c7             	mov    %rax,%rdi
    e3c4:	e8 21 d9 ff ff       	call   bcea <_ZN3MPI4CommD1Ev>
    e3c9:	48 89 d8             	mov    %rbx,%rax
    e3cc:	48 89 c7             	mov    %rax,%rdi
    e3cf:	e8 ac c5 ff ff       	call   a980 <_Unwind_Resume@plt>
    e3d4:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e3d8:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e3df:	00 00 
    e3e1:	74 05                	je     e3e8 <_ZN3MPI9IntracommC1EP19ompi_communicator_t+0xe4>
    e3e3:	e8 e8 c0 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    e3e8:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    e3ec:	c9                   	leave  
    e3ed:	c3                   	ret    

000000000000e3ee <_ZNK3MPI9Intracomm4ScanEPKvPviRKNS_8DatatypeERKNS_2OpE>:
    e3ee:	f3 0f 1e fa          	endbr64 
    e3f2:	55                   	push   %rbp
    e3f3:	48 89 e5             	mov    %rsp,%rbp
    e3f6:	41 54                	push   %r12
    e3f8:	53                   	push   %rbx
    e3f9:	48 83 ec 30          	sub    $0x30,%rsp
    e3fd:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e401:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e405:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    e409:	89 4d d4             	mov    %ecx,-0x2c(%rbp)
    e40c:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    e410:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    e414:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e418:	4c 8b 60 08          	mov    0x8(%rax),%r12
    e41c:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    e420:	48 89 c7             	mov    %rax,%rdi
    e423:	e8 ee 16 00 00       	call   fb16 <_ZNK3MPI2OpcvP9ompi_op_tEv>
    e428:	48 89 c3             	mov    %rax,%rbx
    e42b:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e42f:	48 89 c7             	mov    %rax,%rdi
    e432:	e8 69 d4 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    e437:	48 89 c1             	mov    %rax,%rcx
    e43a:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    e43d:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    e441:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e445:	4d 89 e1             	mov    %r12,%r9
    e448:	49 89 d8             	mov    %rbx,%r8
    e44b:	48 89 c7             	mov    %rax,%rdi
    e44e:	e8 6d c0 ff ff       	call   a4c0 <MPI_Scan@plt>
    e453:	90                   	nop
    e454:	48 83 c4 30          	add    $0x30,%rsp
    e458:	5b                   	pop    %rbx
    e459:	41 5c                	pop    %r12
    e45b:	5d                   	pop    %rbp
    e45c:	c3                   	ret    
    e45d:	90                   	nop

000000000000e45e <_ZNK3MPI9Intracomm6ExscanEPKvPviRKNS_8DatatypeERKNS_2OpE>:
    e45e:	f3 0f 1e fa          	endbr64 
    e462:	55                   	push   %rbp
    e463:	48 89 e5             	mov    %rsp,%rbp
    e466:	41 54                	push   %r12
    e468:	53                   	push   %rbx
    e469:	48 83 ec 30          	sub    $0x30,%rsp
    e46d:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e471:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e475:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    e479:	89 4d d4             	mov    %ecx,-0x2c(%rbp)
    e47c:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    e480:	4c 89 4d c0          	mov    %r9,-0x40(%rbp)
    e484:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e488:	4c 8b 60 08          	mov    0x8(%rax),%r12
    e48c:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    e490:	48 89 c7             	mov    %rax,%rdi
    e493:	e8 7e 16 00 00       	call   fb16 <_ZNK3MPI2OpcvP9ompi_op_tEv>
    e498:	48 89 c3             	mov    %rax,%rbx
    e49b:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e49f:	48 89 c7             	mov    %rax,%rdi
    e4a2:	e8 f9 d3 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    e4a7:	48 89 c1             	mov    %rax,%rcx
    e4aa:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    e4ad:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    e4b1:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e4b5:	4d 89 e1             	mov    %r12,%r9
    e4b8:	49 89 d8             	mov    %rbx,%r8
    e4bb:	48 89 c7             	mov    %rax,%rdi
    e4be:	e8 cd be ff ff       	call   a390 <MPI_Exscan@plt>
    e4c3:	90                   	nop
    e4c4:	48 83 c4 30          	add    $0x30,%rsp
    e4c8:	5b                   	pop    %rbx
    e4c9:	41 5c                	pop    %r12
    e4cb:	5d                   	pop    %rbp
    e4cc:	c3                   	ret    
    e4cd:	90                   	nop

000000000000e4ce <_ZNK3MPI9Intracomm5CloneEv>:
    e4ce:	f3 0f 1e fa          	endbr64 
    e4d2:	55                   	push   %rbp
    e4d3:	48 89 e5             	mov    %rsp,%rbp
    e4d6:	41 54                	push   %r12
    e4d8:	53                   	push   %rbx
    e4d9:	48 83 ec 30          	sub    $0x30,%rsp
    e4dd:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
    e4e1:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e4e8:	00 00 
    e4ea:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    e4ee:	31 c0                	xor    %eax,%eax
    e4f0:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e4f4:	48 8b 40 08          	mov    0x8(%rax),%rax
    e4f8:	48 8d 55 d8          	lea    -0x28(%rbp),%rdx
    e4fc:	48 89 d6             	mov    %rdx,%rsi
    e4ff:	48 89 c7             	mov    %rax,%rdi
    e502:	e8 39 c1 ff ff       	call   a640 <MPI_Comm_dup@plt>
    e507:	bf 10 00 00 00       	mov    $0x10,%edi
    e50c:	e8 bf be ff ff       	call   a3d0 <_Znwm@plt>
    e511:	48 89 c3             	mov    %rax,%rbx
    e514:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e518:	48 89 c6             	mov    %rax,%rsi
    e51b:	48 89 df             	mov    %rbx,%rdi
    e51e:	e8 e1 fd ff ff       	call   e304 <_ZN3MPI9IntracommC1EP19ompi_communicator_t>
    e523:	48 89 5d e0          	mov    %rbx,-0x20(%rbp)
    e527:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e52b:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    e52f:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    e536:	00 00 
    e538:	74 26                	je     e560 <_ZNK3MPI9Intracomm5CloneEv+0x92>
    e53a:	eb 1f                	jmp    e55b <_ZNK3MPI9Intracomm5CloneEv+0x8d>
    e53c:	f3 0f 1e fa          	endbr64 
    e540:	49 89 c4             	mov    %rax,%r12
    e543:	be 10 00 00 00       	mov    $0x10,%esi
    e548:	48 89 df             	mov    %rbx,%rdi
    e54b:	e8 a0 be ff ff       	call   a3f0 <_ZdlPvm@plt>
    e550:	4c 89 e0             	mov    %r12,%rax
    e553:	48 89 c7             	mov    %rax,%rdi
    e556:	e8 25 c4 ff ff       	call   a980 <_Unwind_Resume@plt>
    e55b:	e8 70 bf ff ff       	call   a4d0 <__stack_chk_fail@plt>
    e560:	48 83 c4 30          	add    $0x30,%rsp
    e564:	5b                   	pop    %rbx
    e565:	41 5c                	pop    %r12
    e567:	5d                   	pop    %rbp
    e568:	c3                   	ret    
    e569:	90                   	nop

000000000000e56a <_ZNK3MPI9Intracomm6CreateERKNS_5GroupE>:
    e56a:	f3 0f 1e fa          	endbr64 
    e56e:	55                   	push   %rbp
    e56f:	48 89 e5             	mov    %rsp,%rbp
    e572:	48 83 ec 30          	sub    $0x30,%rsp
    e576:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e57a:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e57e:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    e582:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e589:	00 00 
    e58b:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e58f:	31 c0                	xor    %eax,%eax
    e591:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e595:	48 89 c7             	mov    %rax,%rdi
    e598:	e8 8d d5 ff ff       	call   bb2a <_ZNK3MPI5GroupcvP12ompi_group_tEv>
    e59d:	48 89 c1             	mov    %rax,%rcx
    e5a0:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e5a4:	48 8b 40 08          	mov    0x8(%rax),%rax
    e5a8:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    e5ac:	48 89 ce             	mov    %rcx,%rsi
    e5af:	48 89 c7             	mov    %rax,%rdi
    e5b2:	e8 49 bd ff ff       	call   a300 <MPI_Comm_create@plt>
    e5b7:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    e5bb:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e5bf:	48 89 d6             	mov    %rdx,%rsi
    e5c2:	48 89 c7             	mov    %rax,%rdi
    e5c5:	e8 3a fd ff ff       	call   e304 <_ZN3MPI9IntracommC1EP19ompi_communicator_t>
    e5ca:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e5ce:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e5d5:	00 00 
    e5d7:	74 05                	je     e5de <_ZNK3MPI9Intracomm6CreateERKNS_5GroupE+0x74>
    e5d9:	e8 f2 be ff ff       	call   a4d0 <__stack_chk_fail@plt>
    e5de:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e5e2:	c9                   	leave  
    e5e3:	c3                   	ret    

000000000000e5e4 <_ZNK3MPI9Intracomm5SplitEii>:
    e5e4:	f3 0f 1e fa          	endbr64 
    e5e8:	55                   	push   %rbp
    e5e9:	48 89 e5             	mov    %rsp,%rbp
    e5ec:	48 83 ec 30          	sub    $0x30,%rsp
    e5f0:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e5f4:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e5f8:	89 55 dc             	mov    %edx,-0x24(%rbp)
    e5fb:	89 4d d8             	mov    %ecx,-0x28(%rbp)
    e5fe:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e605:	00 00 
    e607:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e60b:	31 c0                	xor    %eax,%eax
    e60d:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e611:	48 8b 40 08          	mov    0x8(%rax),%rax
    e615:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    e619:	8b 55 d8             	mov    -0x28(%rbp),%edx
    e61c:	8b 75 dc             	mov    -0x24(%rbp),%esi
    e61f:	48 89 c7             	mov    %rax,%rdi
    e622:	e8 39 c2 ff ff       	call   a860 <MPI_Comm_split@plt>
    e627:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    e62b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e62f:	48 89 d6             	mov    %rdx,%rsi
    e632:	48 89 c7             	mov    %rax,%rdi
    e635:	e8 ca fc ff ff       	call   e304 <_ZN3MPI9IntracommC1EP19ompi_communicator_t>
    e63a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e63e:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e645:	00 00 
    e647:	74 05                	je     e64e <_ZNK3MPI9Intracomm5SplitEii+0x6a>
    e649:	e8 82 be ff ff       	call   a4d0 <__stack_chk_fail@plt>
    e64e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e652:	c9                   	leave  
    e653:	c3                   	ret    

000000000000e654 <_ZNK3MPI9Intracomm16Create_intercommEiRKNS_4CommEii>:
    e654:	f3 0f 1e fa          	endbr64 
    e658:	55                   	push   %rbp
    e659:	48 89 e5             	mov    %rsp,%rbp
    e65c:	48 83 ec 40          	sub    $0x40,%rsp
    e660:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e664:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e668:	89 55 dc             	mov    %edx,-0x24(%rbp)
    e66b:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    e66f:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
    e673:	44 89 4d cc          	mov    %r9d,-0x34(%rbp)
    e677:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e67e:	00 00 
    e680:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e684:	31 c0                	xor    %eax,%eax
    e686:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    e68a:	48 89 c7             	mov    %rax,%rdi
    e68d:	e8 2a d5 ff ff       	call   bbbc <_ZNK3MPI9Comm_NullcvP19ompi_communicator_tEv>
    e692:	48 89 c7             	mov    %rax,%rdi
    e695:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e699:	48 8b 40 08          	mov    0x8(%rax),%rax
    e69d:	4c 8d 45 f0          	lea    -0x10(%rbp),%r8
    e6a1:	8b 4d cc             	mov    -0x34(%rbp),%ecx
    e6a4:	8b 55 d8             	mov    -0x28(%rbp),%edx
    e6a7:	8b 75 dc             	mov    -0x24(%rbp),%esi
    e6aa:	4d 89 c1             	mov    %r8,%r9
    e6ad:	41 89 c8             	mov    %ecx,%r8d
    e6b0:	89 d1                	mov    %edx,%ecx
    e6b2:	48 89 fa             	mov    %rdi,%rdx
    e6b5:	48 89 c7             	mov    %rax,%rdi
    e6b8:	e8 63 bd ff ff       	call   a420 <MPI_Intercomm_create@plt>
    e6bd:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    e6c1:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e6c5:	48 89 d6             	mov    %rdx,%rsi
    e6c8:	48 89 c7             	mov    %rax,%rdi
    e6cb:	e8 04 d7 ff ff       	call   bdd4 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
    e6d0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e6d4:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e6db:	00 00 
    e6dd:	74 05                	je     e6e4 <_ZNK3MPI9Intracomm16Create_intercommEiRKNS_4CommEii+0x90>
    e6df:	e8 ec bd ff ff       	call   a4d0 <__stack_chk_fail@plt>
    e6e4:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e6e8:	c9                   	leave  
    e6e9:	c3                   	ret    

000000000000e6ea <_ZNK3MPI9Intracomm11Create_cartEiPKiPKbb>:
    e6ea:	f3 0f 1e fa          	endbr64 
    e6ee:	55                   	push   %rbp
    e6ef:	48 89 e5             	mov    %rsp,%rbp
    e6f2:	48 83 ec 50          	sub    $0x50,%rsp
    e6f6:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    e6fa:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    e6fe:	89 55 cc             	mov    %edx,-0x34(%rbp)
    e701:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    e705:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    e709:	44 89 c8             	mov    %r9d,%eax
    e70c:	88 45 c8             	mov    %al,-0x38(%rbp)
    e70f:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e716:	00 00 
    e718:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e71c:	31 c0                	xor    %eax,%eax
    e71e:	8b 45 cc             	mov    -0x34(%rbp),%eax
    e721:	48 98                	cltq   
    e723:	48 ba fe ff ff ff ff 	movabs $0x1ffffffffffffffe,%rdx
    e72a:	ff ff 1f 
    e72d:	48 39 d0             	cmp    %rdx,%rax
    e730:	77 19                	ja     e74b <_ZNK3MPI9Intracomm11Create_cartEiPKiPKbb+0x61>
    e732:	48 c1 e0 02          	shl    $0x2,%rax
    e736:	48 89 c7             	mov    %rax,%rdi
    e739:	e8 62 b6 ff ff       	call   9da0 <_Znam@plt>
    e73e:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
    e742:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
    e749:	eb 33                	jmp    e77e <_ZNK3MPI9Intracomm11Create_cartEiPKiPKbb+0x94>
    e74b:	e8 90 bf ff ff       	call   a6e0 <__cxa_throw_bad_array_new_length@plt>
    e750:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    e753:	48 63 d0             	movslq %eax,%rdx
    e756:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
    e75a:	48 01 d0             	add    %rdx,%rax
    e75d:	0f b6 00             	movzbl (%rax),%eax
    e760:	8b 55 e4             	mov    -0x1c(%rbp),%edx
    e763:	48 63 d2             	movslq %edx,%rdx
    e766:	48 8d 0c 95 00 00 00 	lea    0x0(,%rdx,4),%rcx
    e76d:	00 
    e76e:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    e772:	48 01 ca             	add    %rcx,%rdx
    e775:	0f b6 c0             	movzbl %al,%eax
    e778:	89 02                	mov    %eax,(%rdx)
    e77a:	83 45 e4 01          	addl   $0x1,-0x1c(%rbp)
    e77e:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    e781:	3b 45 cc             	cmp    -0x34(%rbp),%eax
    e784:	7c ca                	jl     e750 <_ZNK3MPI9Intracomm11Create_cartEiPKiPKbb+0x66>
    e786:	0f b6 7d c8          	movzbl -0x38(%rbp),%edi
    e78a:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    e78e:	48 8b 40 08          	mov    0x8(%rax),%rax
    e792:	4c 8d 45 e8          	lea    -0x18(%rbp),%r8
    e796:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    e79a:	48 8b 55 c0          	mov    -0x40(%rbp),%rdx
    e79e:	8b 75 cc             	mov    -0x34(%rbp),%esi
    e7a1:	4d 89 c1             	mov    %r8,%r9
    e7a4:	41 89 f8             	mov    %edi,%r8d
    e7a7:	48 89 c7             	mov    %rax,%rdi
    e7aa:	e8 61 b6 ff ff       	call   9e10 <MPI_Cart_create@plt>
    e7af:	48 83 7d f0 00       	cmpq   $0x0,-0x10(%rbp)
    e7b4:	74 0c                	je     e7c2 <_ZNK3MPI9Intracomm11Create_cartEiPKiPKbb+0xd8>
    e7b6:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    e7ba:	48 89 c7             	mov    %rax,%rdi
    e7bd:	e8 fe bd ff ff       	call   a5c0 <_ZdaPv@plt>
    e7c2:	48 8d 55 e8          	lea    -0x18(%rbp),%rdx
    e7c6:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e7ca:	48 89 d6             	mov    %rdx,%rsi
    e7cd:	48 89 c7             	mov    %rax,%rdi
    e7d0:	e8 47 05 00 00       	call   ed1c <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t>
    e7d5:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e7d9:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e7e0:	00 00 
    e7e2:	74 05                	je     e7e9 <_ZNK3MPI9Intracomm11Create_cartEiPKiPKbb+0xff>
    e7e4:	e8 e7 bc ff ff       	call   a4d0 <__stack_chk_fail@plt>
    e7e9:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e7ed:	c9                   	leave  
    e7ee:	c3                   	ret    
    e7ef:	90                   	nop

000000000000e7f0 <_ZNK3MPI9Intracomm12Create_graphEiPKiS2_b>:
    e7f0:	f3 0f 1e fa          	endbr64 
    e7f4:	55                   	push   %rbp
    e7f5:	48 89 e5             	mov    %rsp,%rbp
    e7f8:	48 83 ec 40          	sub    $0x40,%rsp
    e7fc:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    e800:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    e804:	89 55 dc             	mov    %edx,-0x24(%rbp)
    e807:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    e80b:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    e80f:	44 89 c8             	mov    %r9d,%eax
    e812:	88 45 d8             	mov    %al,-0x28(%rbp)
    e815:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e81c:	00 00 
    e81e:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    e822:	31 c0                	xor    %eax,%eax
    e824:	0f b6 7d d8          	movzbl -0x28(%rbp),%edi
    e828:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    e82c:	48 8b 40 08          	mov    0x8(%rax),%rax
    e830:	4c 8d 45 f0          	lea    -0x10(%rbp),%r8
    e834:	48 8b 4d c8          	mov    -0x38(%rbp),%rcx
    e838:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    e83c:	8b 75 dc             	mov    -0x24(%rbp),%esi
    e83f:	4d 89 c1             	mov    %r8,%r9
    e842:	41 89 f8             	mov    %edi,%r8d
    e845:	48 89 c7             	mov    %rax,%rdi
    e848:	e8 83 c2 ff ff       	call   aad0 <MPI_Graph_create@plt>
    e84d:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    e851:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e855:	48 89 d6             	mov    %rdx,%rsi
    e858:	48 89 c7             	mov    %rax,%rdi
    e85b:	e8 54 0a 00 00       	call   f2b4 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t>
    e860:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    e864:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e86b:	00 00 
    e86d:	74 05                	je     e874 <_ZNK3MPI9Intracomm12Create_graphEiPKiS2_b+0x84>
    e86f:	e8 5c bc ff ff       	call   a4d0 <__stack_chk_fail@plt>
    e874:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e878:	c9                   	leave  
    e879:	c3                   	ret    

000000000000e87a <_ZNK3MPI9Intracomm6AcceptEPKcRKNS_4InfoEi>:
    e87a:	f3 0f 1e fa          	endbr64 
    e87e:	55                   	push   %rbp
    e87f:	48 89 e5             	mov    %rsp,%rbp
    e882:	53                   	push   %rbx
    e883:	48 83 ec 48          	sub    $0x48,%rsp
    e887:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    e88b:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    e88f:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    e893:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    e897:	44 89 45 bc          	mov    %r8d,-0x44(%rbp)
    e89b:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e8a2:	00 00 
    e8a4:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    e8a8:	31 c0                	xor    %eax,%eax
    e8aa:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    e8ae:	48 8b 58 08          	mov    0x8(%rax),%rbx
    e8b2:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    e8b6:	48 89 c7             	mov    %rax,%rdi
    e8b9:	e8 ca d5 ff ff       	call   be88 <_ZNK3MPI4InfocvP11ompi_info_tEv>
    e8be:	48 89 c6             	mov    %rax,%rsi
    e8c1:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    e8c5:	8b 55 bc             	mov    -0x44(%rbp),%edx
    e8c8:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e8cc:	49 89 c8             	mov    %rcx,%r8
    e8cf:	48 89 d9             	mov    %rbx,%rcx
    e8d2:	48 89 c7             	mov    %rax,%rdi
    e8d5:	e8 d6 c0 ff ff       	call   a9b0 <MPI_Comm_accept@plt>
    e8da:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    e8de:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e8e2:	48 89 d6             	mov    %rdx,%rsi
    e8e5:	48 89 c7             	mov    %rax,%rdi
    e8e8:	e8 e7 d4 ff ff       	call   bdd4 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
    e8ed:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e8f1:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e8f8:	00 00 
    e8fa:	74 05                	je     e901 <_ZNK3MPI9Intracomm6AcceptEPKcRKNS_4InfoEi+0x87>
    e8fc:	e8 cf bb ff ff       	call   a4d0 <__stack_chk_fail@plt>
    e901:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e905:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    e909:	c9                   	leave  
    e90a:	c3                   	ret    
    e90b:	90                   	nop

000000000000e90c <_ZNK3MPI9Intracomm7ConnectEPKcRKNS_4InfoEi>:
    e90c:	f3 0f 1e fa          	endbr64 
    e910:	55                   	push   %rbp
    e911:	48 89 e5             	mov    %rsp,%rbp
    e914:	53                   	push   %rbx
    e915:	48 83 ec 48          	sub    $0x48,%rsp
    e919:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    e91d:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    e921:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    e925:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    e929:	44 89 45 bc          	mov    %r8d,-0x44(%rbp)
    e92d:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e934:	00 00 
    e936:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    e93a:	31 c0                	xor    %eax,%eax
    e93c:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    e940:	48 8b 58 08          	mov    0x8(%rax),%rbx
    e944:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    e948:	48 89 c7             	mov    %rax,%rdi
    e94b:	e8 38 d5 ff ff       	call   be88 <_ZNK3MPI4InfocvP11ompi_info_tEv>
    e950:	48 89 c6             	mov    %rax,%rsi
    e953:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    e957:	8b 55 bc             	mov    -0x44(%rbp),%edx
    e95a:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e95e:	49 89 c8             	mov    %rcx,%r8
    e961:	48 89 d9             	mov    %rbx,%rcx
    e964:	48 89 c7             	mov    %rax,%rdi
    e967:	e8 b4 be ff ff       	call   a820 <MPI_Comm_connect@plt>
    e96c:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    e970:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e974:	48 89 d6             	mov    %rdx,%rsi
    e977:	48 89 c7             	mov    %rax,%rdi
    e97a:	e8 55 d4 ff ff       	call   bdd4 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
    e97f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    e983:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    e98a:	00 00 
    e98c:	74 05                	je     e993 <_ZNK3MPI9Intracomm7ConnectEPKcRKNS_4InfoEi+0x87>
    e98e:	e8 3d bb ff ff       	call   a4d0 <__stack_chk_fail@plt>
    e993:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    e997:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    e99b:	c9                   	leave  
    e99c:	c3                   	ret    
    e99d:	90                   	nop

000000000000e99e <_ZNK3MPI9Intracomm5SpawnEPKcPS2_iRKNS_4InfoEi>:
    e99e:	f3 0f 1e fa          	endbr64 
    e9a2:	55                   	push   %rbp
    e9a3:	48 89 e5             	mov    %rsp,%rbp
    e9a6:	53                   	push   %rbx
    e9a7:	48 83 ec 48          	sub    $0x48,%rsp
    e9ab:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    e9af:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    e9b3:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    e9b7:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    e9bb:	44 89 45 bc          	mov    %r8d,-0x44(%rbp)
    e9bf:	4c 89 4d b0          	mov    %r9,-0x50(%rbp)
    e9c3:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    e9ca:	00 00 
    e9cc:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    e9d0:	31 c0                	xor    %eax,%eax
    e9d2:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    e9d6:	48 8b 58 08          	mov    0x8(%rax),%rbx
    e9da:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
    e9de:	48 89 c7             	mov    %rax,%rdi
    e9e1:	e8 a2 d4 ff ff       	call   be88 <_ZNK3MPI4InfocvP11ompi_info_tEv>
    e9e6:	48 89 c7             	mov    %rax,%rdi
    e9e9:	8b 55 bc             	mov    -0x44(%rbp),%edx
    e9ec:	48 8b 75 c0          	mov    -0x40(%rbp),%rsi
    e9f0:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    e9f4:	6a 00                	push   $0x0
    e9f6:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    e9fa:	51                   	push   %rcx
    e9fb:	49 89 d9             	mov    %rbx,%r9
    e9fe:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    ea02:	48 89 f9             	mov    %rdi,%rcx
    ea05:	48 89 c7             	mov    %rax,%rdi
    ea08:	e8 93 b7 ff ff       	call   a1a0 <MPI_Comm_spawn@plt>
    ea0d:	48 83 c4 10          	add    $0x10,%rsp
    ea11:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    ea15:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ea19:	48 89 d6             	mov    %rdx,%rsi
    ea1c:	48 89 c7             	mov    %rax,%rdi
    ea1f:	e8 b0 d3 ff ff       	call   bdd4 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
    ea24:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ea28:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    ea2f:	00 00 
    ea31:	74 05                	je     ea38 <_ZNK3MPI9Intracomm5SpawnEPKcPS2_iRKNS_4InfoEi+0x9a>
    ea33:	e8 98 ba ff ff       	call   a4d0 <__stack_chk_fail@plt>
    ea38:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ea3c:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    ea40:	c9                   	leave  
    ea41:	c3                   	ret    

000000000000ea42 <_ZNK3MPI9Intracomm5SpawnEPKcPS2_iRKNS_4InfoEiPi>:
    ea42:	f3 0f 1e fa          	endbr64 
    ea46:	55                   	push   %rbp
    ea47:	48 89 e5             	mov    %rsp,%rbp
    ea4a:	53                   	push   %rbx
    ea4b:	48 83 ec 58          	sub    $0x58,%rsp
    ea4f:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    ea53:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    ea57:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    ea5b:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    ea5f:	44 89 45 bc          	mov    %r8d,-0x44(%rbp)
    ea63:	4c 89 4d b0          	mov    %r9,-0x50(%rbp)
    ea67:	48 8b 45 18          	mov    0x18(%rbp),%rax
    ea6b:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
    ea6f:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ea76:	00 00 
    ea78:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    ea7c:	31 c0                	xor    %eax,%eax
    ea7e:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ea82:	48 8b 58 08          	mov    0x8(%rax),%rbx
    ea86:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
    ea8a:	48 89 c7             	mov    %rax,%rdi
    ea8d:	e8 f6 d3 ff ff       	call   be88 <_ZNK3MPI4InfocvP11ompi_info_tEv>
    ea92:	48 89 c7             	mov    %rax,%rdi
    ea95:	8b 55 bc             	mov    -0x44(%rbp),%edx
    ea98:	48 8b 75 c0          	mov    -0x40(%rbp),%rsi
    ea9c:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    eaa0:	ff 75 a8             	push   -0x58(%rbp)
    eaa3:	48 8d 4d e0          	lea    -0x20(%rbp),%rcx
    eaa7:	51                   	push   %rcx
    eaa8:	49 89 d9             	mov    %rbx,%r9
    eaab:	44 8b 45 10          	mov    0x10(%rbp),%r8d
    eaaf:	48 89 f9             	mov    %rdi,%rcx
    eab2:	48 89 c7             	mov    %rax,%rdi
    eab5:	e8 e6 b6 ff ff       	call   a1a0 <MPI_Comm_spawn@plt>
    eaba:	48 83 c4 10          	add    $0x10,%rsp
    eabe:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
    eac2:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    eac6:	48 89 d6             	mov    %rdx,%rsi
    eac9:	48 89 c7             	mov    %rax,%rdi
    eacc:	e8 03 d3 ff ff       	call   bdd4 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
    ead1:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ead5:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    eadc:	00 00 
    eade:	74 05                	je     eae5 <_ZNK3MPI9Intracomm5SpawnEPKcPS2_iRKNS_4InfoEiPi+0xa3>
    eae0:	e8 eb b9 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    eae5:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    eae9:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    eaed:	c9                   	leave  
    eaee:	c3                   	ret    
    eaef:	90                   	nop

000000000000eaf0 <_ZN3MPI9Intracomm14Spawn_multipleEiPPKcPS3_PKiPKNS_4InfoEi>:
    eaf0:	f3 0f 1e fa          	endbr64 
    eaf4:	55                   	push   %rbp
    eaf5:	48 89 e5             	mov    %rsp,%rbp
    eaf8:	48 83 ec 60          	sub    $0x60,%rsp
    eafc:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    eb00:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    eb04:	89 55 cc             	mov    %edx,-0x34(%rbp)
    eb07:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    eb0b:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    eb0f:	4c 89 4d b0          	mov    %r9,-0x50(%rbp)
    eb13:	48 8b 45 10          	mov    0x10(%rbp),%rax
    eb17:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
    eb1b:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    eb22:	00 00 
    eb24:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    eb28:	31 c0                	xor    %eax,%eax
    eb2a:	48 8b 55 a8          	mov    -0x58(%rbp),%rdx
    eb2e:	8b 45 cc             	mov    -0x34(%rbp),%eax
    eb31:	48 89 d6             	mov    %rdx,%rsi
    eb34:	89 c7                	mov    %eax,%edi
    eb36:	e8 7f 00 00 00       	call   ebba <_ZN3MPI9Intracomm24convert_info_to_mpi_infoEiPKNS_4InfoE>
    eb3b:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
    eb3f:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    eb43:	48 8b 78 08          	mov    0x8(%rax),%rdi
    eb47:	44 8b 4d 18          	mov    0x18(%rbp),%r9d
    eb4b:	4c 8b 55 f0          	mov    -0x10(%rbp),%r10
    eb4f:	48 8b 4d b0          	mov    -0x50(%rbp),%rcx
    eb53:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
    eb57:	48 8b 75 c0          	mov    -0x40(%rbp),%rsi
    eb5b:	8b 45 cc             	mov    -0x34(%rbp),%eax
    eb5e:	48 83 ec 08          	sub    $0x8,%rsp
    eb62:	6a 00                	push   $0x0
    eb64:	4c 8d 45 e8          	lea    -0x18(%rbp),%r8
    eb68:	41 50                	push   %r8
    eb6a:	57                   	push   %rdi
    eb6b:	4d 89 d0             	mov    %r10,%r8
    eb6e:	89 c7                	mov    %eax,%edi
    eb70:	e8 1b bd ff ff       	call   a890 <MPI_Comm_spawn_multiple@plt>
    eb75:	48 83 c4 20          	add    $0x20,%rsp
    eb79:	48 83 7d f0 00       	cmpq   $0x0,-0x10(%rbp)
    eb7e:	74 0c                	je     eb8c <_ZN3MPI9Intracomm14Spawn_multipleEiPPKcPS3_PKiPKNS_4InfoEi+0x9c>
    eb80:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    eb84:	48 89 c7             	mov    %rax,%rdi
    eb87:	e8 34 ba ff ff       	call   a5c0 <_ZdaPv@plt>
    eb8c:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    eb90:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    eb94:	48 89 d6             	mov    %rdx,%rsi
    eb97:	48 89 c7             	mov    %rax,%rdi
    eb9a:	e8 35 d2 ff ff       	call   bdd4 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
    eb9f:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    eba3:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    ebaa:	00 00 
    ebac:	74 05                	je     ebb3 <_ZN3MPI9Intracomm14Spawn_multipleEiPPKcPS3_PKiPKNS_4InfoEi+0xc3>
    ebae:	e8 1d b9 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    ebb3:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ebb7:	c9                   	leave  
    ebb8:	c3                   	ret    
    ebb9:	90                   	nop

000000000000ebba <_ZN3MPI9Intracomm24convert_info_to_mpi_infoEiPKNS_4InfoE>:
    ebba:	f3 0f 1e fa          	endbr64 
    ebbe:	55                   	push   %rbp
    ebbf:	48 89 e5             	mov    %rsp,%rbp
    ebc2:	53                   	push   %rbx
    ebc3:	48 83 ec 28          	sub    $0x28,%rsp
    ebc7:	89 7d dc             	mov    %edi,-0x24(%rbp)
    ebca:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    ebce:	8b 45 dc             	mov    -0x24(%rbp),%eax
    ebd1:	48 98                	cltq   
    ebd3:	48 ba ff ff ff ff ff 	movabs $0xfffffffffffffff,%rdx
    ebda:	ff ff 0f 
    ebdd:	48 39 d0             	cmp    %rdx,%rax
    ebe0:	77 19                	ja     ebfb <_ZN3MPI9Intracomm24convert_info_to_mpi_infoEiPKNS_4InfoE+0x41>
    ebe2:	48 c1 e0 03          	shl    $0x3,%rax
    ebe6:	48 89 c7             	mov    %rax,%rdi
    ebe9:	e8 b2 b1 ff ff       	call   9da0 <_Znam@plt>
    ebee:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    ebf2:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
    ebf9:	eb 3d                	jmp    ec38 <_ZN3MPI9Intracomm24convert_info_to_mpi_infoEiPKNS_4InfoE+0x7e>
    ebfb:	e8 e0 ba ff ff       	call   a6e0 <__cxa_throw_bad_array_new_length@plt>
    ec00:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    ec03:	48 98                	cltq   
    ec05:	48 c1 e0 04          	shl    $0x4,%rax
    ec09:	48 89 c2             	mov    %rax,%rdx
    ec0c:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ec10:	48 01 d0             	add    %rdx,%rax
    ec13:	8b 55 e4             	mov    -0x1c(%rbp),%edx
    ec16:	48 63 d2             	movslq %edx,%rdx
    ec19:	48 8d 0c d5 00 00 00 	lea    0x0(,%rdx,8),%rcx
    ec20:	00 
    ec21:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    ec25:	48 8d 1c 11          	lea    (%rcx,%rdx,1),%rbx
    ec29:	48 89 c7             	mov    %rax,%rdi
    ec2c:	e8 57 d2 ff ff       	call   be88 <_ZNK3MPI4InfocvP11ompi_info_tEv>
    ec31:	48 89 03             	mov    %rax,(%rbx)
    ec34:	83 45 e4 01          	addl   $0x1,-0x1c(%rbp)
    ec38:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    ec3b:	3b 45 dc             	cmp    -0x24(%rbp),%eax
    ec3e:	7c c0                	jl     ec00 <_ZN3MPI9Intracomm24convert_info_to_mpi_infoEiPKNS_4InfoE+0x46>
    ec40:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ec44:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    ec48:	c9                   	leave  
    ec49:	c3                   	ret    

000000000000ec4a <_ZN3MPI9Intracomm14Spawn_multipleEiPPKcPS3_PKiPKNS_4InfoEiPi>:
    ec4a:	f3 0f 1e fa          	endbr64 
    ec4e:	55                   	push   %rbp
    ec4f:	48 89 e5             	mov    %rsp,%rbp
    ec52:	48 83 ec 60          	sub    $0x60,%rsp
    ec56:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    ec5a:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    ec5e:	89 55 cc             	mov    %edx,-0x34(%rbp)
    ec61:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    ec65:	4c 89 45 b8          	mov    %r8,-0x48(%rbp)
    ec69:	4c 89 4d b0          	mov    %r9,-0x50(%rbp)
    ec6d:	48 8b 45 10          	mov    0x10(%rbp),%rax
    ec71:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
    ec75:	48 8b 45 20          	mov    0x20(%rbp),%rax
    ec79:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
    ec7d:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ec84:	00 00 
    ec86:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    ec8a:	31 c0                	xor    %eax,%eax
    ec8c:	48 8b 55 a8          	mov    -0x58(%rbp),%rdx
    ec90:	8b 45 cc             	mov    -0x34(%rbp),%eax
    ec93:	48 89 d6             	mov    %rdx,%rsi
    ec96:	89 c7                	mov    %eax,%edi
    ec98:	e8 1d ff ff ff       	call   ebba <_ZN3MPI9Intracomm24convert_info_to_mpi_infoEiPKNS_4InfoE>
    ec9d:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
    eca1:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    eca5:	48 8b 78 08          	mov    0x8(%rax),%rdi
    eca9:	44 8b 4d 18          	mov    0x18(%rbp),%r9d
    ecad:	4c 8b 55 f0          	mov    -0x10(%rbp),%r10
    ecb1:	48 8b 4d b0          	mov    -0x50(%rbp),%rcx
    ecb5:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
    ecb9:	48 8b 75 c0          	mov    -0x40(%rbp),%rsi
    ecbd:	8b 45 cc             	mov    -0x34(%rbp),%eax
    ecc0:	48 83 ec 08          	sub    $0x8,%rsp
    ecc4:	ff 75 a0             	push   -0x60(%rbp)
    ecc7:	4c 8d 45 e8          	lea    -0x18(%rbp),%r8
    eccb:	41 50                	push   %r8
    eccd:	57                   	push   %rdi
    ecce:	4d 89 d0             	mov    %r10,%r8
    ecd1:	89 c7                	mov    %eax,%edi
    ecd3:	e8 b8 bb ff ff       	call   a890 <MPI_Comm_spawn_multiple@plt>
    ecd8:	48 83 c4 20          	add    $0x20,%rsp
    ecdc:	48 83 7d f0 00       	cmpq   $0x0,-0x10(%rbp)
    ece1:	74 0c                	je     ecef <_ZN3MPI9Intracomm14Spawn_multipleEiPPKcPS3_PKiPKNS_4InfoEiPi+0xa5>
    ece3:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    ece7:	48 89 c7             	mov    %rax,%rdi
    ecea:	e8 d1 b8 ff ff       	call   a5c0 <_ZdaPv@plt>
    ecef:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    ecf3:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ecf7:	48 89 d6             	mov    %rdx,%rsi
    ecfa:	48 89 c7             	mov    %rax,%rdi
    ecfd:	e8 d2 d0 ff ff       	call   bdd4 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
    ed02:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ed06:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    ed0d:	00 00 
    ed0f:	74 05                	je     ed16 <_ZN3MPI9Intracomm14Spawn_multipleEiPPKcPS3_PKiPKNS_4InfoEiPi+0xcc>
    ed11:	e8 ba b7 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    ed16:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ed1a:	c9                   	leave  
    ed1b:	c3                   	ret    

000000000000ed1c <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t>:
    ed1c:	f3 0f 1e fa          	endbr64 
    ed20:	55                   	push   %rbp
    ed21:	48 89 e5             	mov    %rsp,%rbp
    ed24:	53                   	push   %rbx
    ed25:	48 83 ec 28          	sub    $0x28,%rsp
    ed29:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    ed2d:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    ed31:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ed38:	00 00 
    ed3a:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    ed3e:	31 c0                	xor    %eax,%eax
    ed40:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ed44:	48 89 c7             	mov    %rax,%rdi
    ed47:	e8 fc cf ff ff       	call   bd48 <_ZN3MPI9IntracommC1Ev>
    ed4c:	48 8d 15 a5 ce 00 00 	lea    0xcea5(%rip),%rdx        # 1bbf8 <_ZTVN3MPI8CartcommE+0x10>
    ed53:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ed57:	48 89 10             	mov    %rdx,(%rax)
    ed5a:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
    ed61:	e8 88 db ff ff       	call   c8ee <_ZN3MPI14Is_initializedEv>
    ed66:	84 c0                	test   %al,%al
    ed68:	74 1a                	je     ed84 <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0x68>
    ed6a:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ed6e:	48 8b 00             	mov    (%rax),%rax
    ed71:	48 8d 15 a8 f0 00 00 	lea    0xf0a8(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    ed78:	48 39 d0             	cmp    %rdx,%rax
    ed7b:	74 07                	je     ed84 <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0x68>
    ed7d:	b8 01 00 00 00       	mov    $0x1,%eax
    ed82:	eb 05                	jmp    ed89 <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0x6d>
    ed84:	b8 00 00 00 00       	mov    $0x0,%eax
    ed89:	84 c0                	test   %al,%al
    ed8b:	74 40                	je     edcd <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0xb1>
    ed8d:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    ed91:	48 8b 00             	mov    (%rax),%rax
    ed94:	48 8d 55 e4          	lea    -0x1c(%rbp),%rdx
    ed98:	48 89 d6             	mov    %rdx,%rsi
    ed9b:	48 89 c7             	mov    %rax,%rdi
    ed9e:	e8 8d b9 ff ff       	call   a730 <MPI_Topo_test@plt>
    eda3:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    eda6:	83 f8 01             	cmp    $0x1,%eax
    eda9:	75 11                	jne    edbc <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0xa0>
    edab:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    edaf:	48 8b 10             	mov    (%rax),%rdx
    edb2:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    edb6:	48 89 50 08          	mov    %rdx,0x8(%rax)
    edba:	eb 40                	jmp    edfc <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0xe0>
    edbc:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    edc0:	48 8d 15 59 f0 00 00 	lea    0xf059(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    edc7:	48 89 50 08          	mov    %rdx,0x8(%rax)
    edcb:	eb 2f                	jmp    edfc <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0xe0>
    edcd:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    edd1:	48 8b 10             	mov    (%rax),%rdx
    edd4:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    edd8:	48 89 50 08          	mov    %rdx,0x8(%rax)
    eddc:	eb 1e                	jmp    edfc <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0xe0>
    edde:	f3 0f 1e fa          	endbr64 
    ede2:	48 89 c3             	mov    %rax,%rbx
    ede5:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    ede9:	48 89 c7             	mov    %rax,%rdi
    edec:	e8 85 cf ff ff       	call   bd76 <_ZN3MPI9IntracommD1Ev>
    edf1:	48 89 d8             	mov    %rbx,%rax
    edf4:	48 89 c7             	mov    %rax,%rdi
    edf7:	e8 84 bb ff ff       	call   a980 <_Unwind_Resume@plt>
    edfc:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ee00:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    ee07:	00 00 
    ee09:	74 05                	je     ee10 <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t+0xf4>
    ee0b:	e8 c0 b6 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    ee10:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    ee14:	c9                   	leave  
    ee15:	c3                   	ret    

000000000000ee16 <_ZNK3MPI8Cartcomm7Get_dimEv>:
    ee16:	f3 0f 1e fa          	endbr64 
    ee1a:	55                   	push   %rbp
    ee1b:	48 89 e5             	mov    %rsp,%rbp
    ee1e:	48 83 ec 20          	sub    $0x20,%rsp
    ee22:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    ee26:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ee2d:	00 00 
    ee2f:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    ee33:	31 c0                	xor    %eax,%eax
    ee35:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ee39:	48 8b 40 08          	mov    0x8(%rax),%rax
    ee3d:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    ee41:	48 89 d6             	mov    %rdx,%rsi
    ee44:	48 89 c7             	mov    %rax,%rdi
    ee47:	e8 e4 b5 ff ff       	call   a430 <MPI_Cartdim_get@plt>
    ee4c:	8b 45 f4             	mov    -0xc(%rbp),%eax
    ee4f:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    ee53:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    ee5a:	00 00 
    ee5c:	74 05                	je     ee63 <_ZNK3MPI8Cartcomm7Get_dimEv+0x4d>
    ee5e:	e8 6d b6 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    ee63:	c9                   	leave  
    ee64:	c3                   	ret    
    ee65:	90                   	nop

000000000000ee66 <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_>:
    ee66:	f3 0f 1e fa          	endbr64 
    ee6a:	55                   	push   %rbp
    ee6b:	48 89 e5             	mov    %rsp,%rbp
    ee6e:	48 83 ec 40          	sub    $0x40,%rsp
    ee72:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    ee76:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    ee79:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    ee7d:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    ee81:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    ee85:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    ee88:	48 98                	cltq   
    ee8a:	48 ba fe ff ff ff ff 	movabs $0x1ffffffffffffffe,%rdx
    ee91:	ff ff 1f 
    ee94:	48 39 d0             	cmp    %rdx,%rax
    ee97:	77 19                	ja     eeb2 <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_+0x4c>
    ee99:	48 c1 e0 02          	shl    $0x2,%rax
    ee9d:	48 89 c7             	mov    %rax,%rdi
    eea0:	e8 fb ae ff ff       	call   9da0 <_Znam@plt>
    eea5:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    eea9:	c7 45 f4 00 00 00 00 	movl   $0x0,-0xc(%rbp)
    eeb0:	eb 33                	jmp    eee5 <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_+0x7f>
    eeb2:	e8 29 b8 ff ff       	call   a6e0 <__cxa_throw_bad_array_new_length@plt>
    eeb7:	8b 45 f4             	mov    -0xc(%rbp),%eax
    eeba:	48 63 d0             	movslq %eax,%rdx
    eebd:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    eec1:	48 01 d0             	add    %rdx,%rax
    eec4:	0f b6 00             	movzbl (%rax),%eax
    eec7:	8b 55 f4             	mov    -0xc(%rbp),%edx
    eeca:	48 63 d2             	movslq %edx,%rdx
    eecd:	48 8d 0c 95 00 00 00 	lea    0x0(,%rdx,4),%rcx
    eed4:	00 
    eed5:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    eed9:	48 01 ca             	add    %rcx,%rdx
    eedc:	0f b6 c0             	movzbl %al,%eax
    eedf:	89 02                	mov    %eax,(%rdx)
    eee1:	83 45 f4 01          	addl   $0x1,-0xc(%rbp)
    eee5:	8b 45 f4             	mov    -0xc(%rbp),%eax
    eee8:	3b 45 e4             	cmp    -0x1c(%rbp),%eax
    eeeb:	7c ca                	jl     eeb7 <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_+0x51>
    eeed:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    eef1:	48 8b 40 08          	mov    0x8(%rax),%rax
    eef5:	48 8b 7d c8          	mov    -0x38(%rbp),%rdi
    eef9:	48 8b 4d f8          	mov    -0x8(%rbp),%rcx
    eefd:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
    ef01:	8b 75 e4             	mov    -0x1c(%rbp),%esi
    ef04:	49 89 f8             	mov    %rdi,%r8
    ef07:	48 89 c7             	mov    %rax,%rdi
    ef0a:	e8 01 b4 ff ff       	call   a310 <MPI_Cart_get@plt>
    ef0f:	c7 45 f4 00 00 00 00 	movl   $0x0,-0xc(%rbp)
    ef16:	eb 2e                	jmp    ef46 <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_+0xe0>
    ef18:	8b 45 f4             	mov    -0xc(%rbp),%eax
    ef1b:	48 98                	cltq   
    ef1d:	48 8d 14 85 00 00 00 	lea    0x0(,%rax,4),%rdx
    ef24:	00 
    ef25:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ef29:	48 01 d0             	add    %rdx,%rax
    ef2c:	8b 00                	mov    (%rax),%eax
    ef2e:	8b 55 f4             	mov    -0xc(%rbp),%edx
    ef31:	48 63 ca             	movslq %edx,%rcx
    ef34:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    ef38:	48 01 ca             	add    %rcx,%rdx
    ef3b:	85 c0                	test   %eax,%eax
    ef3d:	0f 95 c0             	setne  %al
    ef40:	88 02                	mov    %al,(%rdx)
    ef42:	83 45 f4 01          	addl   $0x1,-0xc(%rbp)
    ef46:	8b 45 f4             	mov    -0xc(%rbp),%eax
    ef49:	3b 45 e4             	cmp    -0x1c(%rbp),%eax
    ef4c:	7c ca                	jl     ef18 <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_+0xb2>
    ef4e:	48 83 7d f8 00       	cmpq   $0x0,-0x8(%rbp)
    ef53:	74 0c                	je     ef61 <_ZNK3MPI8Cartcomm8Get_topoEiPiPbS1_+0xfb>
    ef55:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ef59:	48 89 c7             	mov    %rax,%rdi
    ef5c:	e8 5f b6 ff ff       	call   a5c0 <_ZdaPv@plt>
    ef61:	90                   	nop
    ef62:	c9                   	leave  
    ef63:	c3                   	ret    

000000000000ef64 <_ZNK3MPI8Cartcomm13Get_cart_rankEPKi>:
    ef64:	f3 0f 1e fa          	endbr64 
    ef68:	55                   	push   %rbp
    ef69:	48 89 e5             	mov    %rsp,%rbp
    ef6c:	48 83 ec 20          	sub    $0x20,%rsp
    ef70:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    ef74:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    ef78:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ef7f:	00 00 
    ef81:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    ef85:	31 c0                	xor    %eax,%eax
    ef87:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ef8b:	48 8b 40 08          	mov    0x8(%rax),%rax
    ef8f:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    ef93:	48 8b 4d e0          	mov    -0x20(%rbp),%rcx
    ef97:	48 89 ce             	mov    %rcx,%rsi
    ef9a:	48 89 c7             	mov    %rax,%rdi
    ef9d:	e8 fe b9 ff ff       	call   a9a0 <MPI_Cart_rank@plt>
    efa2:	8b 45 f4             	mov    -0xc(%rbp),%eax
    efa5:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    efa9:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    efb0:	00 00 
    efb2:	74 05                	je     efb9 <_ZNK3MPI8Cartcomm13Get_cart_rankEPKi+0x55>
    efb4:	e8 17 b5 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    efb9:	c9                   	leave  
    efba:	c3                   	ret    
    efbb:	90                   	nop

000000000000efbc <_ZNK3MPI8Cartcomm10Get_coordsEiiPi>:
    efbc:	f3 0f 1e fa          	endbr64 
    efc0:	55                   	push   %rbp
    efc1:	48 89 e5             	mov    %rsp,%rbp
    efc4:	48 83 ec 20          	sub    $0x20,%rsp
    efc8:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    efcc:	89 75 f4             	mov    %esi,-0xc(%rbp)
    efcf:	89 55 f0             	mov    %edx,-0x10(%rbp)
    efd2:	48 89 4d e8          	mov    %rcx,-0x18(%rbp)
    efd6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    efda:	48 8b 40 08          	mov    0x8(%rax),%rax
    efde:	48 8b 4d e8          	mov    -0x18(%rbp),%rcx
    efe2:	8b 55 f0             	mov    -0x10(%rbp),%edx
    efe5:	8b 75 f4             	mov    -0xc(%rbp),%esi
    efe8:	48 89 c7             	mov    %rax,%rdi
    efeb:	e8 50 b4 ff ff       	call   a440 <MPI_Cart_coords@plt>
    eff0:	90                   	nop
    eff1:	c9                   	leave  
    eff2:	c3                   	ret    
    eff3:	90                   	nop

000000000000eff4 <_ZNK3MPI8Cartcomm5ShiftEiiRiS1_>:
    eff4:	f3 0f 1e fa          	endbr64 
    eff8:	55                   	push   %rbp
    eff9:	48 89 e5             	mov    %rsp,%rbp
    effc:	48 83 ec 20          	sub    $0x20,%rsp
    f000:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    f004:	89 75 f4             	mov    %esi,-0xc(%rbp)
    f007:	89 55 f0             	mov    %edx,-0x10(%rbp)
    f00a:	48 89 4d e8          	mov    %rcx,-0x18(%rbp)
    f00e:	4c 89 45 e0          	mov    %r8,-0x20(%rbp)
    f012:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f016:	48 8b 40 08          	mov    0x8(%rax),%rax
    f01a:	48 8b 7d e0          	mov    -0x20(%rbp),%rdi
    f01e:	48 8b 4d e8          	mov    -0x18(%rbp),%rcx
    f022:	8b 55 f0             	mov    -0x10(%rbp),%edx
    f025:	8b 75 f4             	mov    -0xc(%rbp),%esi
    f028:	49 89 f8             	mov    %rdi,%r8
    f02b:	48 89 c7             	mov    %rax,%rdi
    f02e:	e8 8d b0 ff ff       	call   a0c0 <MPI_Cart_shift@plt>
    f033:	90                   	nop
    f034:	c9                   	leave  
    f035:	c3                   	ret    

000000000000f036 <_ZNK3MPI8Cartcomm3SubEPKb>:
    f036:	f3 0f 1e fa          	endbr64 
    f03a:	55                   	push   %rbp
    f03b:	48 89 e5             	mov    %rsp,%rbp
    f03e:	48 83 ec 40          	sub    $0x40,%rsp
    f042:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    f046:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    f04a:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    f04e:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f055:	00 00 
    f057:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f05b:	31 c0                	xor    %eax,%eax
    f05d:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    f061:	48 8b 40 08          	mov    0x8(%rax),%rax
    f065:	48 8d 55 e0          	lea    -0x20(%rbp),%rdx
    f069:	48 89 d6             	mov    %rdx,%rsi
    f06c:	48 89 c7             	mov    %rax,%rdi
    f06f:	e8 bc b3 ff ff       	call   a430 <MPI_Cartdim_get@plt>
    f074:	8b 45 e0             	mov    -0x20(%rbp),%eax
    f077:	48 98                	cltq   
    f079:	48 ba fe ff ff ff ff 	movabs $0x1ffffffffffffffe,%rdx
    f080:	ff ff 1f 
    f083:	48 39 d0             	cmp    %rdx,%rax
    f086:	77 19                	ja     f0a1 <_ZNK3MPI8Cartcomm3SubEPKb+0x6b>
    f088:	48 c1 e0 02          	shl    $0x2,%rax
    f08c:	48 89 c7             	mov    %rax,%rdi
    f08f:	e8 0c ad ff ff       	call   9da0 <_Znam@plt>
    f094:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
    f098:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
    f09f:	eb 33                	jmp    f0d4 <_ZNK3MPI8Cartcomm3SubEPKb+0x9e>
    f0a1:	e8 3a b6 ff ff       	call   a6e0 <__cxa_throw_bad_array_new_length@plt>
    f0a6:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    f0a9:	48 63 d0             	movslq %eax,%rdx
    f0ac:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    f0b0:	48 01 d0             	add    %rdx,%rax
    f0b3:	0f b6 00             	movzbl (%rax),%eax
    f0b6:	8b 55 e4             	mov    -0x1c(%rbp),%edx
    f0b9:	48 63 d2             	movslq %edx,%rdx
    f0bc:	48 8d 0c 95 00 00 00 	lea    0x0(,%rdx,4),%rcx
    f0c3:	00 
    f0c4:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f0c8:	48 01 ca             	add    %rcx,%rdx
    f0cb:	0f b6 c0             	movzbl %al,%eax
    f0ce:	89 02                	mov    %eax,(%rdx)
    f0d0:	83 45 e4 01          	addl   $0x1,-0x1c(%rbp)
    f0d4:	8b 45 e0             	mov    -0x20(%rbp),%eax
    f0d7:	39 45 e4             	cmp    %eax,-0x1c(%rbp)
    f0da:	7c ca                	jl     f0a6 <_ZNK3MPI8Cartcomm3SubEPKb+0x70>
    f0dc:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    f0e0:	48 8b 40 08          	mov    0x8(%rax),%rax
    f0e4:	48 8d 55 e8          	lea    -0x18(%rbp),%rdx
    f0e8:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    f0ec:	48 89 ce             	mov    %rcx,%rsi
    f0ef:	48 89 c7             	mov    %rax,%rdi
    f0f2:	e8 29 b0 ff ff       	call   a120 <MPI_Cart_sub@plt>
    f0f7:	48 83 7d f0 00       	cmpq   $0x0,-0x10(%rbp)
    f0fc:	74 0c                	je     f10a <_ZNK3MPI8Cartcomm3SubEPKb+0xd4>
    f0fe:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    f102:	48 89 c7             	mov    %rax,%rdi
    f105:	e8 b6 b4 ff ff       	call   a5c0 <_ZdaPv@plt>
    f10a:	48 8d 55 e8          	lea    -0x18(%rbp),%rdx
    f10e:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f112:	48 89 d6             	mov    %rdx,%rsi
    f115:	48 89 c7             	mov    %rax,%rdi
    f118:	e8 ff fb ff ff       	call   ed1c <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t>
    f11d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f121:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f128:	00 00 
    f12a:	74 05                	je     f131 <_ZNK3MPI8Cartcomm3SubEPKb+0xfb>
    f12c:	e8 9f b3 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f131:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f135:	c9                   	leave  
    f136:	c3                   	ret    
    f137:	90                   	nop

000000000000f138 <_ZNK3MPI8Cartcomm3MapEiPKiPKb>:
    f138:	f3 0f 1e fa          	endbr64 
    f13c:	55                   	push   %rbp
    f13d:	48 89 e5             	mov    %rsp,%rbp
    f140:	48 83 ec 40          	sub    $0x40,%rsp
    f144:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    f148:	89 75 d4             	mov    %esi,-0x2c(%rbp)
    f14b:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
    f14f:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
    f153:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f15a:	00 00 
    f15c:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f160:	31 c0                	xor    %eax,%eax
    f162:	8b 45 d4             	mov    -0x2c(%rbp),%eax
    f165:	48 98                	cltq   
    f167:	48 ba fe ff ff ff ff 	movabs $0x1ffffffffffffffe,%rdx
    f16e:	ff ff 1f 
    f171:	48 39 d0             	cmp    %rdx,%rax
    f174:	77 19                	ja     f18f <_ZNK3MPI8Cartcomm3MapEiPKiPKb+0x57>
    f176:	48 c1 e0 02          	shl    $0x2,%rax
    f17a:	48 89 c7             	mov    %rax,%rdi
    f17d:	e8 1e ac ff ff       	call   9da0 <_Znam@plt>
    f182:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
    f186:	c7 45 ec 00 00 00 00 	movl   $0x0,-0x14(%rbp)
    f18d:	eb 33                	jmp    f1c2 <_ZNK3MPI8Cartcomm3MapEiPKiPKb+0x8a>
    f18f:	e8 4c b5 ff ff       	call   a6e0 <__cxa_throw_bad_array_new_length@plt>
    f194:	8b 45 ec             	mov    -0x14(%rbp),%eax
    f197:	48 63 d0             	movslq %eax,%rdx
    f19a:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
    f19e:	48 01 d0             	add    %rdx,%rax
    f1a1:	0f b6 00             	movzbl (%rax),%eax
    f1a4:	8b 55 ec             	mov    -0x14(%rbp),%edx
    f1a7:	48 63 d2             	movslq %edx,%rdx
    f1aa:	48 8d 0c 95 00 00 00 	lea    0x0(,%rdx,4),%rcx
    f1b1:	00 
    f1b2:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f1b6:	48 01 ca             	add    %rcx,%rdx
    f1b9:	0f b6 c0             	movzbl %al,%eax
    f1bc:	89 02                	mov    %eax,(%rdx)
    f1be:	83 45 ec 01          	addl   $0x1,-0x14(%rbp)
    f1c2:	8b 45 ec             	mov    -0x14(%rbp),%eax
    f1c5:	3b 45 d4             	cmp    -0x2c(%rbp),%eax
    f1c8:	7c ca                	jl     f194 <_ZNK3MPI8Cartcomm3MapEiPKiPKb+0x5c>
    f1ca:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f1ce:	48 8b 40 08          	mov    0x8(%rax),%rax
    f1d2:	48 8d 7d e8          	lea    -0x18(%rbp),%rdi
    f1d6:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    f1da:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
    f1de:	8b 75 d4             	mov    -0x2c(%rbp),%esi
    f1e1:	49 89 f8             	mov    %rdi,%r8
    f1e4:	48 89 c7             	mov    %rax,%rdi
    f1e7:	e8 94 b1 ff ff       	call   a380 <MPI_Cart_map@plt>
    f1ec:	48 83 7d f0 00       	cmpq   $0x0,-0x10(%rbp)
    f1f1:	74 0c                	je     f1ff <_ZNK3MPI8Cartcomm3MapEiPKiPKb+0xc7>
    f1f3:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    f1f7:	48 89 c7             	mov    %rax,%rdi
    f1fa:	e8 c1 b3 ff ff       	call   a5c0 <_ZdaPv@plt>
    f1ff:	8b 45 e8             	mov    -0x18(%rbp),%eax
    f202:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    f206:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f20d:	00 00 
    f20f:	74 05                	je     f216 <_ZNK3MPI8Cartcomm3MapEiPKiPKb+0xde>
    f211:	e8 ba b2 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f216:	c9                   	leave  
    f217:	c3                   	ret    

000000000000f218 <_ZNK3MPI8Cartcomm5CloneEv>:
    f218:	f3 0f 1e fa          	endbr64 
    f21c:	55                   	push   %rbp
    f21d:	48 89 e5             	mov    %rsp,%rbp
    f220:	41 54                	push   %r12
    f222:	53                   	push   %rbx
    f223:	48 83 ec 30          	sub    $0x30,%rsp
    f227:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
    f22b:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f232:	00 00 
    f234:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    f238:	31 c0                	xor    %eax,%eax
    f23a:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    f23e:	48 8b 40 08          	mov    0x8(%rax),%rax
    f242:	48 8d 55 d8          	lea    -0x28(%rbp),%rdx
    f246:	48 89 d6             	mov    %rdx,%rsi
    f249:	48 89 c7             	mov    %rax,%rdi
    f24c:	e8 ef b3 ff ff       	call   a640 <MPI_Comm_dup@plt>
    f251:	bf 10 00 00 00       	mov    $0x10,%edi
    f256:	e8 75 b1 ff ff       	call   a3d0 <_Znwm@plt>
    f25b:	48 89 c3             	mov    %rax,%rbx
    f25e:	48 8d 45 d8          	lea    -0x28(%rbp),%rax
    f262:	48 89 c6             	mov    %rax,%rsi
    f265:	48 89 df             	mov    %rbx,%rdi
    f268:	e8 af fa ff ff       	call   ed1c <_ZN3MPI8CartcommC1ERKP19ompi_communicator_t>
    f26d:	48 89 5d e0          	mov    %rbx,-0x20(%rbp)
    f271:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f275:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    f279:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f280:	00 00 
    f282:	74 26                	je     f2aa <_ZNK3MPI8Cartcomm5CloneEv+0x92>
    f284:	eb 1f                	jmp    f2a5 <_ZNK3MPI8Cartcomm5CloneEv+0x8d>
    f286:	f3 0f 1e fa          	endbr64 
    f28a:	49 89 c4             	mov    %rax,%r12
    f28d:	be 10 00 00 00       	mov    $0x10,%esi
    f292:	48 89 df             	mov    %rbx,%rdi
    f295:	e8 56 b1 ff ff       	call   a3f0 <_ZdlPvm@plt>
    f29a:	4c 89 e0             	mov    %r12,%rax
    f29d:	48 89 c7             	mov    %rax,%rdi
    f2a0:	e8 db b6 ff ff       	call   a980 <_Unwind_Resume@plt>
    f2a5:	e8 26 b2 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f2aa:	48 83 c4 30          	add    $0x30,%rsp
    f2ae:	5b                   	pop    %rbx
    f2af:	41 5c                	pop    %r12
    f2b1:	5d                   	pop    %rbp
    f2b2:	c3                   	ret    
    f2b3:	90                   	nop

000000000000f2b4 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t>:
    f2b4:	f3 0f 1e fa          	endbr64 
    f2b8:	55                   	push   %rbp
    f2b9:	48 89 e5             	mov    %rsp,%rbp
    f2bc:	53                   	push   %rbx
    f2bd:	48 83 ec 28          	sub    $0x28,%rsp
    f2c1:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
    f2c5:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
    f2c9:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f2d0:	00 00 
    f2d2:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    f2d6:	31 c0                	xor    %eax,%eax
    f2d8:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f2dc:	48 89 c7             	mov    %rax,%rdi
    f2df:	e8 64 ca ff ff       	call   bd48 <_ZN3MPI9IntracommC1Ev>
    f2e4:	48 8d 15 ad c6 00 00 	lea    0xc6ad(%rip),%rdx        # 1b998 <_ZTVN3MPI9GraphcommE+0x10>
    f2eb:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f2ef:	48 89 10             	mov    %rdx,(%rax)
    f2f2:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
    f2f9:	e8 f0 d5 ff ff       	call   c8ee <_ZN3MPI14Is_initializedEv>
    f2fe:	84 c0                	test   %al,%al
    f300:	74 1a                	je     f31c <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0x68>
    f302:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    f306:	48 8b 00             	mov    (%rax),%rax
    f309:	48 8d 15 10 eb 00 00 	lea    0xeb10(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    f310:	48 39 d0             	cmp    %rdx,%rax
    f313:	74 07                	je     f31c <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0x68>
    f315:	b8 01 00 00 00       	mov    $0x1,%eax
    f31a:	eb 05                	jmp    f321 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0x6d>
    f31c:	b8 00 00 00 00       	mov    $0x0,%eax
    f321:	84 c0                	test   %al,%al
    f323:	74 40                	je     f365 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0xb1>
    f325:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    f329:	48 8b 00             	mov    (%rax),%rax
    f32c:	48 8d 55 e4          	lea    -0x1c(%rbp),%rdx
    f330:	48 89 d6             	mov    %rdx,%rsi
    f333:	48 89 c7             	mov    %rax,%rdi
    f336:	e8 f5 b3 ff ff       	call   a730 <MPI_Topo_test@plt>
    f33b:	8b 45 e4             	mov    -0x1c(%rbp),%eax
    f33e:	83 f8 02             	cmp    $0x2,%eax
    f341:	75 11                	jne    f354 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0xa0>
    f343:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    f347:	48 8b 10             	mov    (%rax),%rdx
    f34a:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f34e:	48 89 50 08          	mov    %rdx,0x8(%rax)
    f352:	eb 40                	jmp    f394 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0xe0>
    f354:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f358:	48 8d 15 c1 ea 00 00 	lea    0xeac1(%rip),%rdx        # 1de20 <ompi_mpi_comm_null@@Base>
    f35f:	48 89 50 08          	mov    %rdx,0x8(%rax)
    f363:	eb 2f                	jmp    f394 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0xe0>
    f365:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
    f369:	48 8b 10             	mov    (%rax),%rdx
    f36c:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f370:	48 89 50 08          	mov    %rdx,0x8(%rax)
    f374:	eb 1e                	jmp    f394 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0xe0>
    f376:	f3 0f 1e fa          	endbr64 
    f37a:	48 89 c3             	mov    %rax,%rbx
    f37d:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f381:	48 89 c7             	mov    %rax,%rdi
    f384:	e8 ed c9 ff ff       	call   bd76 <_ZN3MPI9IntracommD1Ev>
    f389:	48 89 d8             	mov    %rbx,%rax
    f38c:	48 89 c7             	mov    %rax,%rdi
    f38f:	e8 ec b5 ff ff       	call   a980 <_Unwind_Resume@plt>
    f394:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f398:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f39f:	00 00 
    f3a1:	74 05                	je     f3a8 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t+0xf4>
    f3a3:	e8 28 b1 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f3a8:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    f3ac:	c9                   	leave  
    f3ad:	c3                   	ret    

000000000000f3ae <_ZNK3MPI9Graphcomm8Get_dimsEPiS1_>:
    f3ae:	f3 0f 1e fa          	endbr64 
    f3b2:	55                   	push   %rbp
    f3b3:	48 89 e5             	mov    %rsp,%rbp
    f3b6:	48 83 ec 20          	sub    $0x20,%rsp
    f3ba:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    f3be:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    f3c2:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
    f3c6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f3ca:	48 8b 40 08          	mov    0x8(%rax),%rax
    f3ce:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    f3d2:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
    f3d6:	48 89 ce             	mov    %rcx,%rsi
    f3d9:	48 89 c7             	mov    %rax,%rdi
    f3dc:	e8 bf ac ff ff       	call   a0a0 <MPI_Graphdims_get@plt>
    f3e1:	90                   	nop
    f3e2:	c9                   	leave  
    f3e3:	c3                   	ret    

000000000000f3e4 <_ZNK3MPI9Graphcomm8Get_topoEiiPiS1_>:
    f3e4:	f3 0f 1e fa          	endbr64 
    f3e8:	55                   	push   %rbp
    f3e9:	48 89 e5             	mov    %rsp,%rbp
    f3ec:	48 83 ec 20          	sub    $0x20,%rsp
    f3f0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    f3f4:	89 75 f4             	mov    %esi,-0xc(%rbp)
    f3f7:	89 55 f0             	mov    %edx,-0x10(%rbp)
    f3fa:	48 89 4d e8          	mov    %rcx,-0x18(%rbp)
    f3fe:	4c 89 45 e0          	mov    %r8,-0x20(%rbp)
    f402:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f406:	48 8b 40 08          	mov    0x8(%rax),%rax
    f40a:	48 8b 7d e0          	mov    -0x20(%rbp),%rdi
    f40e:	48 8b 4d e8          	mov    -0x18(%rbp),%rcx
    f412:	8b 55 f0             	mov    -0x10(%rbp),%edx
    f415:	8b 75 f4             	mov    -0xc(%rbp),%esi
    f418:	49 89 f8             	mov    %rdi,%r8
    f41b:	48 89 c7             	mov    %rax,%rdi
    f41e:	e8 bd ae ff ff       	call   a2e0 <MPI_Graph_get@plt>
    f423:	90                   	nop
    f424:	c9                   	leave  
    f425:	c3                   	ret    

000000000000f426 <_ZNK3MPI9Graphcomm19Get_neighbors_countEi>:
    f426:	f3 0f 1e fa          	endbr64 
    f42a:	55                   	push   %rbp
    f42b:	48 89 e5             	mov    %rsp,%rbp
    f42e:	48 83 ec 20          	sub    $0x20,%rsp
    f432:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f436:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    f439:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f440:	00 00 
    f442:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f446:	31 c0                	xor    %eax,%eax
    f448:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f44c:	48 8b 40 08          	mov    0x8(%rax),%rax
    f450:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    f454:	8b 4d e4             	mov    -0x1c(%rbp),%ecx
    f457:	89 ce                	mov    %ecx,%esi
    f459:	48 89 c7             	mov    %rax,%rdi
    f45c:	e8 0f b6 ff ff       	call   aa70 <MPI_Graph_neighbors_count@plt>
    f461:	8b 45 f4             	mov    -0xc(%rbp),%eax
    f464:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    f468:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f46f:	00 00 
    f471:	74 05                	je     f478 <_ZNK3MPI9Graphcomm19Get_neighbors_countEi+0x52>
    f473:	e8 58 b0 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f478:	c9                   	leave  
    f479:	c3                   	ret    

000000000000f47a <_ZNK3MPI9Graphcomm13Get_neighborsEiiPi>:
    f47a:	f3 0f 1e fa          	endbr64 
    f47e:	55                   	push   %rbp
    f47f:	48 89 e5             	mov    %rsp,%rbp
    f482:	48 83 ec 20          	sub    $0x20,%rsp
    f486:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    f48a:	89 75 f4             	mov    %esi,-0xc(%rbp)
    f48d:	89 55 f0             	mov    %edx,-0x10(%rbp)
    f490:	48 89 4d e8          	mov    %rcx,-0x18(%rbp)
    f494:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f498:	48 8b 40 08          	mov    0x8(%rax),%rax
    f49c:	48 8b 4d e8          	mov    -0x18(%rbp),%rcx
    f4a0:	8b 55 f0             	mov    -0x10(%rbp),%edx
    f4a3:	8b 75 f4             	mov    -0xc(%rbp),%esi
    f4a6:	48 89 c7             	mov    %rax,%rdi
    f4a9:	e8 82 aa ff ff       	call   9f30 <MPI_Graph_neighbors@plt>
    f4ae:	90                   	nop
    f4af:	c9                   	leave  
    f4b0:	c3                   	ret    
    f4b1:	90                   	nop

000000000000f4b2 <_ZNK3MPI9Graphcomm3MapEiPKiS2_>:
    f4b2:	f3 0f 1e fa          	endbr64 
    f4b6:	55                   	push   %rbp
    f4b7:	48 89 e5             	mov    %rsp,%rbp
    f4ba:	48 83 ec 30          	sub    $0x30,%rsp
    f4be:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f4c2:	89 75 e4             	mov    %esi,-0x1c(%rbp)
    f4c5:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    f4c9:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    f4cd:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f4d4:	00 00 
    f4d6:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f4da:	31 c0                	xor    %eax,%eax
    f4dc:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f4e0:	48 8b 40 08          	mov    0x8(%rax),%rax
    f4e4:	48 8d 7d f4          	lea    -0xc(%rbp),%rdi
    f4e8:	48 8b 4d d0          	mov    -0x30(%rbp),%rcx
    f4ec:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
    f4f0:	8b 75 e4             	mov    -0x1c(%rbp),%esi
    f4f3:	49 89 f8             	mov    %rdi,%r8
    f4f6:	48 89 c7             	mov    %rax,%rdi
    f4f9:	e8 52 ae ff ff       	call   a350 <MPI_Graph_map@plt>
    f4fe:	8b 45 f4             	mov    -0xc(%rbp),%eax
    f501:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    f505:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f50c:	00 00 
    f50e:	74 05                	je     f515 <_ZNK3MPI9Graphcomm3MapEiPKiS2_+0x63>
    f510:	e8 bb af ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f515:	c9                   	leave  
    f516:	c3                   	ret    
    f517:	90                   	nop

000000000000f518 <_ZNK3MPI9Graphcomm5CloneEv>:
    f518:	f3 0f 1e fa          	endbr64 
    f51c:	55                   	push   %rbp
    f51d:	48 89 e5             	mov    %rsp,%rbp
    f520:	41 54                	push   %r12
    f522:	53                   	push   %rbx
    f523:	48 83 ec 30          	sub    $0x30,%rsp
    f527:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
    f52b:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f532:	00 00 
    f534:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    f538:	31 c0                	xor    %eax,%eax
    f53a:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    f53e:	48 8b 40 08          	mov    0x8(%rax),%rax
    f542:	48 8d 55 d8          	lea    -0x28(%rbp),%rdx
    f546:	48 89 d6             	mov    %rdx,%rsi
    f549:	48 89 c7             	mov    %rax,%rdi
    f54c:	e8 ef b0 ff ff       	call   a640 <MPI_Comm_dup@plt>
    f551:	bf 10 00 00 00       	mov    $0x10,%edi
    f556:	e8 75 ae ff ff       	call   a3d0 <_Znwm@plt>
    f55b:	48 89 c3             	mov    %rax,%rbx
    f55e:	48 8d 45 d8          	lea    -0x28(%rbp),%rax
    f562:	48 89 c6             	mov    %rax,%rsi
    f565:	48 89 df             	mov    %rbx,%rdi
    f568:	e8 47 fd ff ff       	call   f2b4 <_ZN3MPI9GraphcommC1ERKP19ompi_communicator_t>
    f56d:	48 89 5d e0          	mov    %rbx,-0x20(%rbp)
    f571:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f575:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    f579:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f580:	00 00 
    f582:	74 26                	je     f5aa <_ZNK3MPI9Graphcomm5CloneEv+0x92>
    f584:	eb 1f                	jmp    f5a5 <_ZNK3MPI9Graphcomm5CloneEv+0x8d>
    f586:	f3 0f 1e fa          	endbr64 
    f58a:	49 89 c4             	mov    %rax,%r12
    f58d:	be 10 00 00 00       	mov    $0x10,%esi
    f592:	48 89 df             	mov    %rbx,%rdi
    f595:	e8 56 ae ff ff       	call   a3f0 <_ZdlPvm@plt>
    f59a:	4c 89 e0             	mov    %r12,%rax
    f59d:	48 89 c7             	mov    %rax,%rdi
    f5a0:	e8 db b3 ff ff       	call   a980 <_Unwind_Resume@plt>
    f5a5:	e8 26 af ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f5aa:	48 83 c4 30          	add    $0x30,%rsp
    f5ae:	5b                   	pop    %rbx
    f5af:	41 5c                	pop    %r12
    f5b1:	5d                   	pop    %rbp
    f5b2:	c3                   	ret    
    f5b3:	90                   	nop

000000000000f5b4 <_ZNK3MPI9Intercomm5CloneEv>:
    f5b4:	f3 0f 1e fa          	endbr64 
    f5b8:	55                   	push   %rbp
    f5b9:	48 89 e5             	mov    %rsp,%rbp
    f5bc:	53                   	push   %rbx
    f5bd:	48 83 ec 38          	sub    $0x38,%rsp
    f5c1:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
    f5c5:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f5cc:	00 00 
    f5ce:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
    f5d2:	31 c0                	xor    %eax,%eax
    f5d4:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    f5d8:	48 8b 40 08          	mov    0x8(%rax),%rax
    f5dc:	48 8d 55 d8          	lea    -0x28(%rbp),%rdx
    f5e0:	48 89 d6             	mov    %rdx,%rsi
    f5e3:	48 89 c7             	mov    %rax,%rdi
    f5e6:	e8 55 b0 ff ff       	call   a640 <MPI_Comm_dup@plt>
    f5eb:	bf 10 00 00 00       	mov    $0x10,%edi
    f5f0:	e8 db ad ff ff       	call   a3d0 <_Znwm@plt>
    f5f5:	48 89 c3             	mov    %rax,%rbx
    f5f8:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f5fc:	48 89 c6             	mov    %rax,%rsi
    f5ff:	48 89 df             	mov    %rbx,%rdi
    f602:	e8 cd c7 ff ff       	call   bdd4 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
    f607:	48 89 5d e0          	mov    %rbx,-0x20(%rbp)
    f60b:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f60f:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
    f613:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f61a:	00 00 
    f61c:	74 05                	je     f623 <_ZNK3MPI9Intercomm5CloneEv+0x6f>
    f61e:	e8 ad ae ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f623:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    f627:	c9                   	leave  
    f628:	c3                   	ret    
    f629:	90                   	nop

000000000000f62a <_ZNK3MPI9Intercomm15Get_remote_sizeEv>:
    f62a:	f3 0f 1e fa          	endbr64 
    f62e:	55                   	push   %rbp
    f62f:	48 89 e5             	mov    %rsp,%rbp
    f632:	48 83 ec 20          	sub    $0x20,%rsp
    f636:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f63a:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f641:	00 00 
    f643:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f647:	31 c0                	xor    %eax,%eax
    f649:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f64d:	48 8b 40 08          	mov    0x8(%rax),%rax
    f651:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    f655:	48 89 d6             	mov    %rdx,%rsi
    f658:	48 89 c7             	mov    %rax,%rdi
    f65b:	e8 30 b1 ff ff       	call   a790 <MPI_Comm_remote_size@plt>
    f660:	8b 45 f4             	mov    -0xc(%rbp),%eax
    f663:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    f667:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f66e:	00 00 
    f670:	74 05                	je     f677 <_ZNK3MPI9Intercomm15Get_remote_sizeEv+0x4d>
    f672:	e8 59 ae ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f677:	c9                   	leave  
    f678:	c3                   	ret    
    f679:	90                   	nop

000000000000f67a <_ZNK3MPI9Intercomm16Get_remote_groupEv>:
    f67a:	f3 0f 1e fa          	endbr64 
    f67e:	55                   	push   %rbp
    f67f:	48 89 e5             	mov    %rsp,%rbp
    f682:	48 83 ec 20          	sub    $0x20,%rsp
    f686:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f68a:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f68e:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f695:	00 00 
    f697:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f69b:	31 c0                	xor    %eax,%eax
    f69d:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f6a1:	48 8b 40 08          	mov    0x8(%rax),%rax
    f6a5:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    f6a9:	48 89 d6             	mov    %rdx,%rsi
    f6ac:	48 89 c7             	mov    %rax,%rdi
    f6af:	e8 5c af ff ff       	call   a610 <MPI_Comm_remote_group@plt>
    f6b4:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f6b8:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f6bc:	48 89 d6             	mov    %rdx,%rsi
    f6bf:	48 89 c7             	mov    %rax,%rdi
    f6c2:	e8 e7 c3 ff ff       	call   baae <_ZN3MPI5GroupC1EP12ompi_group_t>
    f6c7:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f6cb:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f6d2:	00 00 
    f6d4:	74 05                	je     f6db <_ZNK3MPI9Intercomm16Get_remote_groupEv+0x61>
    f6d6:	e8 f5 ad ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f6db:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f6df:	c9                   	leave  
    f6e0:	c3                   	ret    
    f6e1:	90                   	nop

000000000000f6e2 <_ZNK3MPI9Intercomm5MergeEb>:
    f6e2:	f3 0f 1e fa          	endbr64 
    f6e6:	55                   	push   %rbp
    f6e7:	48 89 e5             	mov    %rsp,%rbp
    f6ea:	48 83 ec 30          	sub    $0x30,%rsp
    f6ee:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f6f2:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f6f6:	89 d0                	mov    %edx,%eax
    f6f8:	88 45 dc             	mov    %al,-0x24(%rbp)
    f6fb:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f702:	00 00 
    f704:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f708:	31 c0                	xor    %eax,%eax
    f70a:	0f b6 4d dc          	movzbl -0x24(%rbp),%ecx
    f70e:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f712:	48 8b 40 08          	mov    0x8(%rax),%rax
    f716:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    f71a:	89 ce                	mov    %ecx,%esi
    f71c:	48 89 c7             	mov    %rax,%rdi
    f71f:	e8 ec b2 ff ff       	call   aa10 <MPI_Intercomm_merge@plt>
    f724:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f728:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f72c:	48 89 d6             	mov    %rdx,%rsi
    f72f:	48 89 c7             	mov    %rax,%rdi
    f732:	e8 cd eb ff ff       	call   e304 <_ZN3MPI9IntracommC1EP19ompi_communicator_t>
    f737:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f73b:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f742:	00 00 
    f744:	74 05                	je     f74b <_ZNK3MPI9Intercomm5MergeEb+0x69>
    f746:	e8 85 ad ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f74b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f74f:	c9                   	leave  
    f750:	c3                   	ret    
    f751:	90                   	nop

000000000000f752 <_ZNK3MPI9Intercomm6CreateERKNS_5GroupE>:
    f752:	f3 0f 1e fa          	endbr64 
    f756:	55                   	push   %rbp
    f757:	48 89 e5             	mov    %rsp,%rbp
    f75a:	48 83 ec 30          	sub    $0x30,%rsp
    f75e:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f762:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f766:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    f76a:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f771:	00 00 
    f773:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f777:	31 c0                	xor    %eax,%eax
    f779:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
    f77d:	48 89 c7             	mov    %rax,%rdi
    f780:	e8 a5 c3 ff ff       	call   bb2a <_ZNK3MPI5GroupcvP12ompi_group_tEv>
    f785:	48 89 c1             	mov    %rax,%rcx
    f788:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f78c:	48 8b 40 08          	mov    0x8(%rax),%rax
    f790:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    f794:	48 89 ce             	mov    %rcx,%rsi
    f797:	48 89 c7             	mov    %rax,%rdi
    f79a:	e8 61 ab ff ff       	call   a300 <MPI_Comm_create@plt>
    f79f:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f7a3:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f7a7:	48 89 d6             	mov    %rdx,%rsi
    f7aa:	48 89 c7             	mov    %rax,%rdi
    f7ad:	e8 22 c6 ff ff       	call   bdd4 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
    f7b2:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f7b6:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f7bd:	00 00 
    f7bf:	74 05                	je     f7c6 <_ZNK3MPI9Intercomm6CreateERKNS_5GroupE+0x74>
    f7c1:	e8 0a ad ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f7c6:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f7ca:	c9                   	leave  
    f7cb:	c3                   	ret    

000000000000f7cc <_ZNK3MPI9Intercomm5SplitEii>:
    f7cc:	f3 0f 1e fa          	endbr64 
    f7d0:	55                   	push   %rbp
    f7d1:	48 89 e5             	mov    %rsp,%rbp
    f7d4:	48 83 ec 30          	sub    $0x30,%rsp
    f7d8:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f7dc:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f7e0:	89 55 dc             	mov    %edx,-0x24(%rbp)
    f7e3:	89 4d d8             	mov    %ecx,-0x28(%rbp)
    f7e6:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f7ed:	00 00 
    f7ef:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f7f3:	31 c0                	xor    %eax,%eax
    f7f5:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f7f9:	48 8b 40 08          	mov    0x8(%rax),%rax
    f7fd:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    f801:	8b 55 d8             	mov    -0x28(%rbp),%edx
    f804:	8b 75 dc             	mov    -0x24(%rbp),%esi
    f807:	48 89 c7             	mov    %rax,%rdi
    f80a:	e8 51 b0 ff ff       	call   a860 <MPI_Comm_split@plt>
    f80f:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f813:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f817:	48 89 d6             	mov    %rdx,%rsi
    f81a:	48 89 c7             	mov    %rax,%rdi
    f81d:	e8 b2 c5 ff ff       	call   bdd4 <_ZN3MPI9IntercommC1EP19ompi_communicator_t>
    f822:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f826:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f82d:	00 00 
    f82f:	74 05                	je     f836 <_ZNK3MPI9Intercomm5SplitEii+0x6a>
    f831:	e8 9a ac ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f836:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f83a:	c9                   	leave  
    f83b:	c3                   	ret    

000000000000f83c <_ZNK3MPI5Group8Get_sizeEv>:
    f83c:	f3 0f 1e fa          	endbr64 
    f840:	55                   	push   %rbp
    f841:	48 89 e5             	mov    %rsp,%rbp
    f844:	48 83 ec 20          	sub    $0x20,%rsp
    f848:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f84c:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f853:	00 00 
    f855:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f859:	31 c0                	xor    %eax,%eax
    f85b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f85f:	48 8b 40 08          	mov    0x8(%rax),%rax
    f863:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    f867:	48 89 d6             	mov    %rdx,%rsi
    f86a:	48 89 c7             	mov    %rax,%rdi
    f86d:	e8 3e a6 ff ff       	call   9eb0 <MPI_Group_size@plt>
    f872:	8b 45 f4             	mov    -0xc(%rbp),%eax
    f875:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    f879:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f880:	00 00 
    f882:	74 05                	je     f889 <_ZNK3MPI5Group8Get_sizeEv+0x4d>
    f884:	e8 47 ac ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f889:	c9                   	leave  
    f88a:	c3                   	ret    
    f88b:	90                   	nop

000000000000f88c <_ZNK3MPI5Group8Get_rankEv>:
    f88c:	f3 0f 1e fa          	endbr64 
    f890:	55                   	push   %rbp
    f891:	48 89 e5             	mov    %rsp,%rbp
    f894:	48 83 ec 20          	sub    $0x20,%rsp
    f898:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f89c:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f8a3:	00 00 
    f8a5:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f8a9:	31 c0                	xor    %eax,%eax
    f8ab:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f8af:	48 8b 40 08          	mov    0x8(%rax),%rax
    f8b3:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    f8b7:	48 89 d6             	mov    %rdx,%rsi
    f8ba:	48 89 c7             	mov    %rax,%rdi
    f8bd:	e8 3e b1 ff ff       	call   aa00 <MPI_Group_rank@plt>
    f8c2:	8b 45 f4             	mov    -0xc(%rbp),%eax
    f8c5:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    f8c9:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    f8d0:	00 00 
    f8d2:	74 05                	je     f8d9 <_ZNK3MPI5Group8Get_rankEv+0x4d>
    f8d4:	e8 f7 ab ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f8d9:	c9                   	leave  
    f8da:	c3                   	ret    
    f8db:	90                   	nop

000000000000f8dc <_ZNK3MPI5Group4InclEiPKi>:
    f8dc:	f3 0f 1e fa          	endbr64 
    f8e0:	55                   	push   %rbp
    f8e1:	48 89 e5             	mov    %rsp,%rbp
    f8e4:	48 83 ec 30          	sub    $0x30,%rsp
    f8e8:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f8ec:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f8f0:	89 55 dc             	mov    %edx,-0x24(%rbp)
    f8f3:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    f8f7:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f8fe:	00 00 
    f900:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f904:	31 c0                	xor    %eax,%eax
    f906:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f90a:	48 8b 40 08          	mov    0x8(%rax),%rax
    f90e:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    f912:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    f916:	8b 75 dc             	mov    -0x24(%rbp),%esi
    f919:	48 89 c7             	mov    %rax,%rdi
    f91c:	e8 ef a8 ff ff       	call   a210 <MPI_Group_incl@plt>
    f921:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f925:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f929:	48 89 d6             	mov    %rdx,%rsi
    f92c:	48 89 c7             	mov    %rax,%rdi
    f92f:	e8 7a c1 ff ff       	call   baae <_ZN3MPI5GroupC1EP12ompi_group_t>
    f934:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f938:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f93f:	00 00 
    f941:	74 05                	je     f948 <_ZNK3MPI5Group4InclEiPKi+0x6c>
    f943:	e8 88 ab ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f948:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f94c:	c9                   	leave  
    f94d:	c3                   	ret    

000000000000f94e <_ZNK3MPI5Group4ExclEiPKi>:
    f94e:	f3 0f 1e fa          	endbr64 
    f952:	55                   	push   %rbp
    f953:	48 89 e5             	mov    %rsp,%rbp
    f956:	48 83 ec 30          	sub    $0x30,%rsp
    f95a:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f95e:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f962:	89 55 dc             	mov    %edx,-0x24(%rbp)
    f965:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    f969:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f970:	00 00 
    f972:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f976:	31 c0                	xor    %eax,%eax
    f978:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f97c:	48 8b 40 08          	mov    0x8(%rax),%rax
    f980:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    f984:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    f988:	8b 75 dc             	mov    -0x24(%rbp),%esi
    f98b:	48 89 c7             	mov    %rax,%rdi
    f98e:	e8 4d a4 ff ff       	call   9de0 <MPI_Group_excl@plt>
    f993:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    f997:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f99b:	48 89 d6             	mov    %rdx,%rsi
    f99e:	48 89 c7             	mov    %rax,%rdi
    f9a1:	e8 08 c1 ff ff       	call   baae <_ZN3MPI5GroupC1EP12ompi_group_t>
    f9a6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    f9aa:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    f9b1:	00 00 
    f9b3:	74 05                	je     f9ba <_ZNK3MPI5Group4ExclEiPKi+0x6c>
    f9b5:	e8 16 ab ff ff       	call   a4d0 <__stack_chk_fail@plt>
    f9ba:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    f9be:	c9                   	leave  
    f9bf:	c3                   	ret    

000000000000f9c0 <_ZNK3MPI5Group10Range_inclEiPA3_Ki>:
    f9c0:	f3 0f 1e fa          	endbr64 
    f9c4:	55                   	push   %rbp
    f9c5:	48 89 e5             	mov    %rsp,%rbp
    f9c8:	48 83 ec 30          	sub    $0x30,%rsp
    f9cc:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    f9d0:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    f9d4:	89 55 dc             	mov    %edx,-0x24(%rbp)
    f9d7:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    f9db:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    f9e2:	00 00 
    f9e4:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    f9e8:	31 c0                	xor    %eax,%eax
    f9ea:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    f9ee:	48 8b 40 08          	mov    0x8(%rax),%rax
    f9f2:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    f9f6:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    f9fa:	8b 75 dc             	mov    -0x24(%rbp),%esi
    f9fd:	48 89 c7             	mov    %rax,%rdi
    fa00:	e8 7b a4 ff ff       	call   9e80 <MPI_Group_range_incl@plt>
    fa05:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    fa09:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fa0d:	48 89 d6             	mov    %rdx,%rsi
    fa10:	48 89 c7             	mov    %rax,%rdi
    fa13:	e8 96 c0 ff ff       	call   baae <_ZN3MPI5GroupC1EP12ompi_group_t>
    fa18:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fa1c:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    fa23:	00 00 
    fa25:	74 05                	je     fa2c <_ZNK3MPI5Group10Range_inclEiPA3_Ki+0x6c>
    fa27:	e8 a4 aa ff ff       	call   a4d0 <__stack_chk_fail@plt>
    fa2c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fa30:	c9                   	leave  
    fa31:	c3                   	ret    

000000000000fa32 <_ZNK3MPI5Group10Range_exclEiPA3_Ki>:
    fa32:	f3 0f 1e fa          	endbr64 
    fa36:	55                   	push   %rbp
    fa37:	48 89 e5             	mov    %rsp,%rbp
    fa3a:	48 83 ec 30          	sub    $0x30,%rsp
    fa3e:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fa42:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    fa46:	89 55 dc             	mov    %edx,-0x24(%rbp)
    fa49:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    fa4d:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    fa54:	00 00 
    fa56:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    fa5a:	31 c0                	xor    %eax,%eax
    fa5c:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    fa60:	48 8b 40 08          	mov    0x8(%rax),%rax
    fa64:	48 8d 4d f0          	lea    -0x10(%rbp),%rcx
    fa68:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
    fa6c:	8b 75 dc             	mov    -0x24(%rbp),%esi
    fa6f:	48 89 c7             	mov    %rax,%rdi
    fa72:	e8 09 ab ff ff       	call   a580 <MPI_Group_range_excl@plt>
    fa77:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    fa7b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fa7f:	48 89 d6             	mov    %rdx,%rsi
    fa82:	48 89 c7             	mov    %rax,%rdi
    fa85:	e8 24 c0 ff ff       	call   baae <_ZN3MPI5GroupC1EP12ompi_group_t>
    fa8a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fa8e:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    fa95:	00 00 
    fa97:	74 05                	je     fa9e <_ZNK3MPI5Group10Range_exclEiPA3_Ki+0x6c>
    fa99:	e8 32 aa ff ff       	call   a4d0 <__stack_chk_fail@plt>
    fa9e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    faa2:	c9                   	leave  
    faa3:	c3                   	ret    

000000000000faa4 <_ZN3MPI5Group4FreeEv>:
    faa4:	f3 0f 1e fa          	endbr64 
    faa8:	55                   	push   %rbp
    faa9:	48 89 e5             	mov    %rsp,%rbp
    faac:	48 83 ec 10          	sub    $0x10,%rsp
    fab0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fab4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fab8:	48 83 c0 08          	add    $0x8,%rax
    fabc:	48 89 c7             	mov    %rax,%rdi
    fabf:	e8 9c ae ff ff       	call   a960 <MPI_Group_free@plt>
    fac4:	90                   	nop
    fac5:	c9                   	leave  
    fac6:	c3                   	ret    
    fac7:	90                   	nop

000000000000fac8 <_ZN3MPI2OpD1Ev>:
    fac8:	f3 0f 1e fa          	endbr64 
    facc:	55                   	push   %rbp
    facd:	48 89 e5             	mov    %rsp,%rbp
    fad0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fad4:	48 8d 15 2d bc 00 00 	lea    0xbc2d(%rip),%rdx        # 1b708 <_ZTVN3MPI2OpE+0x10>
    fadb:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fadf:	48 89 10             	mov    %rdx,(%rax)
    fae2:	90                   	nop
    fae3:	5d                   	pop    %rbp
    fae4:	c3                   	ret    
    fae5:	90                   	nop

000000000000fae6 <_ZN3MPI2OpD0Ev>:
    fae6:	f3 0f 1e fa          	endbr64 
    faea:	55                   	push   %rbp
    faeb:	48 89 e5             	mov    %rsp,%rbp
    faee:	48 83 ec 10          	sub    $0x10,%rsp
    faf2:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    faf6:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fafa:	48 89 c7             	mov    %rax,%rdi
    fafd:	e8 c6 ff ff ff       	call   fac8 <_ZN3MPI2OpD1Ev>
    fb02:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fb06:	be 10 00 00 00       	mov    $0x10,%esi
    fb0b:	48 89 c7             	mov    %rax,%rdi
    fb0e:	e8 dd a8 ff ff       	call   a3f0 <_ZdlPvm@plt>
    fb13:	c9                   	leave  
    fb14:	c3                   	ret    
    fb15:	90                   	nop

000000000000fb16 <_ZNK3MPI2OpcvP9ompi_op_tEv>:
    fb16:	f3 0f 1e fa          	endbr64 
    fb1a:	55                   	push   %rbp
    fb1b:	48 89 e5             	mov    %rsp,%rbp
    fb1e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fb22:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fb26:	48 8b 40 08          	mov    0x8(%rax),%rax
    fb2a:	5d                   	pop    %rbp
    fb2b:	c3                   	ret    

000000000000fb2c <_ZN3MPI2Op4InitEPFvPKvPviRKNS_8DatatypeEEb>:
    fb2c:	f3 0f 1e fa          	endbr64 
    fb30:	55                   	push   %rbp
    fb31:	48 89 e5             	mov    %rsp,%rbp
    fb34:	48 83 ec 20          	sub    $0x20,%rsp
    fb38:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fb3c:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    fb40:	89 d0                	mov    %edx,%eax
    fb42:	88 45 ec             	mov    %al,-0x14(%rbp)
    fb45:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fb49:	48 8d 50 08          	lea    0x8(%rax),%rdx
    fb4d:	0f b6 45 ec          	movzbl -0x14(%rbp),%eax
    fb51:	89 c6                	mov    %eax,%esi
    fb53:	48 8b 05 66 d4 00 00 	mov    0xd466(%rip),%rax        # 1cfc0 <ompi_mpi_cxx_op_intercept@Base>
    fb5a:	48 89 c7             	mov    %rax,%rdi
    fb5d:	e8 3e a3 ff ff       	call   9ea0 <MPI_Op_create@plt>
    fb62:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fb66:	48 8b 40 08          	mov    0x8(%rax),%rax
    fb6a:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    fb6e:	48 89 d6             	mov    %rdx,%rsi
    fb71:	48 89 c7             	mov    %rax,%rdi
    fb74:	e8 d7 a8 ff ff       	call   a450 <ompi_op_set_cxx_callback@plt>
    fb79:	90                   	nop
    fb7a:	c9                   	leave  
    fb7b:	c3                   	ret    

000000000000fb7c <_ZN3MPI2Op4FreeEv>:
    fb7c:	f3 0f 1e fa          	endbr64 
    fb80:	55                   	push   %rbp
    fb81:	48 89 e5             	mov    %rsp,%rbp
    fb84:	48 83 ec 10          	sub    $0x10,%rsp
    fb88:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fb8c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fb90:	48 83 c0 08          	add    $0x8,%rax
    fb94:	48 89 c7             	mov    %rax,%rdi
    fb97:	e8 64 ac ff ff       	call   a800 <MPI_Op_free@plt>
    fb9c:	90                   	nop
    fb9d:	c9                   	leave  
    fb9e:	c3                   	ret    
    fb9f:	90                   	nop

000000000000fba0 <_ZNK3MPI2Op12Reduce_localEPKvPviRKNS_8DatatypeE>:
    fba0:	f3 0f 1e fa          	endbr64 
    fba4:	55                   	push   %rbp
    fba5:	48 89 e5             	mov    %rsp,%rbp
    fba8:	53                   	push   %rbx
    fba9:	48 83 ec 38          	sub    $0x38,%rsp
    fbad:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fbb1:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    fbb5:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
    fbb9:	89 4d d4             	mov    %ecx,-0x2c(%rbp)
    fbbc:	4c 89 45 c8          	mov    %r8,-0x38(%rbp)
    fbc0:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fbc4:	48 8b 58 08          	mov    0x8(%rax),%rbx
    fbc8:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
    fbcc:	48 89 c7             	mov    %rax,%rdi
    fbcf:	e8 cc bc ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    fbd4:	48 89 c1             	mov    %rax,%rcx
    fbd7:	8b 55 d4             	mov    -0x2c(%rbp),%edx
    fbda:	48 8b 75 d8          	mov    -0x28(%rbp),%rsi
    fbde:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    fbe2:	49 89 d8             	mov    %rbx,%r8
    fbe5:	48 89 c7             	mov    %rax,%rdi
    fbe8:	e8 23 a9 ff ff       	call   a510 <MPI_Reduce_local@plt>
    fbed:	90                   	nop
    fbee:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
    fbf2:	c9                   	leave  
    fbf3:	c3                   	ret    

000000000000fbf4 <_ZNK3MPI2Op14Is_commutativeEv>:
    fbf4:	f3 0f 1e fa          	endbr64 
    fbf8:	55                   	push   %rbp
    fbf9:	48 89 e5             	mov    %rsp,%rbp
    fbfc:	48 83 ec 20          	sub    $0x20,%rsp
    fc00:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fc04:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    fc0b:	00 00 
    fc0d:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    fc11:	31 c0                	xor    %eax,%eax
    fc13:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fc17:	48 8b 40 08          	mov    0x8(%rax),%rax
    fc1b:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    fc1f:	48 89 d6             	mov    %rdx,%rsi
    fc22:	48 89 c7             	mov    %rax,%rdi
    fc25:	e8 16 a7 ff ff       	call   a340 <MPI_Op_commutative@plt>
    fc2a:	8b 45 f4             	mov    -0xc(%rbp),%eax
    fc2d:	85 c0                	test   %eax,%eax
    fc2f:	0f 95 c0             	setne  %al
    fc32:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    fc36:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    fc3d:	00 00 
    fc3f:	74 05                	je     fc46 <_ZNK3MPI2Op14Is_commutativeEv+0x52>
    fc41:	e8 8a a8 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    fc46:	c9                   	leave  
    fc47:	c3                   	ret    

000000000000fc48 <_ZN3MPI10Errhandler4FreeEv>:
    fc48:	f3 0f 1e fa          	endbr64 
    fc4c:	55                   	push   %rbp
    fc4d:	48 89 e5             	mov    %rsp,%rbp
    fc50:	48 83 ec 10          	sub    $0x10,%rsp
    fc54:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fc58:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fc5c:	48 83 c0 08          	add    $0x8,%rax
    fc60:	48 89 c7             	mov    %rax,%rdi
    fc63:	e8 78 a7 ff ff       	call   a3e0 <MPI_Errhandler_free@plt>
    fc68:	90                   	nop
    fc69:	c9                   	leave  
    fc6a:	c3                   	ret    
    fc6b:	90                   	nop

000000000000fc6c <_ZNK3MPI6Status9Get_countERKNS_8DatatypeE>:
    fc6c:	f3 0f 1e fa          	endbr64 
    fc70:	55                   	push   %rbp
    fc71:	48 89 e5             	mov    %rsp,%rbp
    fc74:	48 83 ec 20          	sub    $0x20,%rsp
    fc78:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fc7c:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    fc80:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    fc87:	00 00 
    fc89:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    fc8d:	31 c0                	xor    %eax,%eax
    fc8f:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    fc93:	48 89 c7             	mov    %rax,%rdi
    fc96:	e8 05 bc ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    fc9b:	48 89 c6             	mov    %rax,%rsi
    fc9e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fca2:	48 8d 48 08          	lea    0x8(%rax),%rcx
    fca6:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
    fcaa:	48 89 c2             	mov    %rax,%rdx
    fcad:	48 89 cf             	mov    %rcx,%rdi
    fcb0:	e8 8b aa ff ff       	call   a740 <MPI_Get_count@plt>
    fcb5:	8b 45 f4             	mov    -0xc(%rbp),%eax
    fcb8:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    fcbc:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    fcc3:	00 00 
    fcc5:	74 05                	je     fccc <_ZNK3MPI6Status9Get_countERKNS_8DatatypeE+0x60>
    fcc7:	e8 04 a8 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    fccc:	c9                   	leave  
    fccd:	c3                   	ret    

000000000000fcce <_ZNK3MPI6Status12Is_cancelledEv>:
    fcce:	f3 0f 1e fa          	endbr64 
    fcd2:	55                   	push   %rbp
    fcd3:	48 89 e5             	mov    %rsp,%rbp
    fcd6:	48 83 ec 20          	sub    $0x20,%rsp
    fcda:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fcde:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    fce5:	00 00 
    fce7:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    fceb:	31 c0                	xor    %eax,%eax
    fced:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fcf1:	48 8d 50 08          	lea    0x8(%rax),%rdx
    fcf5:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
    fcf9:	48 89 c6             	mov    %rax,%rsi
    fcfc:	48 89 d7             	mov    %rdx,%rdi
    fcff:	e8 6c a6 ff ff       	call   a370 <MPI_Test_cancelled@plt>
    fd04:	8b 45 f4             	mov    -0xc(%rbp),%eax
    fd07:	85 c0                	test   %eax,%eax
    fd09:	0f 95 c0             	setne  %al
    fd0c:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    fd10:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    fd17:	00 00 
    fd19:	74 05                	je     fd20 <_ZNK3MPI6Status12Is_cancelledEv+0x52>
    fd1b:	e8 b0 a7 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    fd20:	c9                   	leave  
    fd21:	c3                   	ret    

000000000000fd22 <_ZNK3MPI6Status12Get_elementsERKNS_8DatatypeE>:
    fd22:	f3 0f 1e fa          	endbr64 
    fd26:	55                   	push   %rbp
    fd27:	48 89 e5             	mov    %rsp,%rbp
    fd2a:	48 83 ec 20          	sub    $0x20,%rsp
    fd2e:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fd32:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    fd36:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    fd3d:	00 00 
    fd3f:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    fd43:	31 c0                	xor    %eax,%eax
    fd45:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    fd49:	48 89 c7             	mov    %rax,%rdi
    fd4c:	e8 4f bb ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    fd51:	48 89 c6             	mov    %rax,%rsi
    fd54:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fd58:	48 8d 48 08          	lea    0x8(%rax),%rcx
    fd5c:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
    fd60:	48 89 c2             	mov    %rax,%rdx
    fd63:	48 89 cf             	mov    %rcx,%rdi
    fd66:	e8 05 a7 ff ff       	call   a470 <MPI_Get_elements@plt>
    fd6b:	8b 45 f4             	mov    -0xc(%rbp),%eax
    fd6e:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    fd72:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    fd79:	00 00 
    fd7b:	74 05                	je     fd82 <_ZNK3MPI6Status12Get_elementsERKNS_8DatatypeE+0x60>
    fd7d:	e8 4e a7 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    fd82:	c9                   	leave  
    fd83:	c3                   	ret    

000000000000fd84 <_ZNK3MPI6Status10Get_sourceEv>:
    fd84:	f3 0f 1e fa          	endbr64 
    fd88:	55                   	push   %rbp
    fd89:	48 89 e5             	mov    %rsp,%rbp
    fd8c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fd90:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fd94:	8b 40 08             	mov    0x8(%rax),%eax
    fd97:	89 45 fc             	mov    %eax,-0x4(%rbp)
    fd9a:	8b 45 fc             	mov    -0x4(%rbp),%eax
    fd9d:	5d                   	pop    %rbp
    fd9e:	c3                   	ret    
    fd9f:	90                   	nop

000000000000fda0 <_ZN3MPI6Status10Set_sourceEi>:
    fda0:	f3 0f 1e fa          	endbr64 
    fda4:	55                   	push   %rbp
    fda5:	48 89 e5             	mov    %rsp,%rbp
    fda8:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fdac:	89 75 f4             	mov    %esi,-0xc(%rbp)
    fdaf:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fdb3:	8b 55 f4             	mov    -0xc(%rbp),%edx
    fdb6:	89 50 08             	mov    %edx,0x8(%rax)
    fdb9:	90                   	nop
    fdba:	5d                   	pop    %rbp
    fdbb:	c3                   	ret    

000000000000fdbc <_ZNK3MPI6Status7Get_tagEv>:
    fdbc:	f3 0f 1e fa          	endbr64 
    fdc0:	55                   	push   %rbp
    fdc1:	48 89 e5             	mov    %rsp,%rbp
    fdc4:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fdc8:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fdcc:	8b 40 0c             	mov    0xc(%rax),%eax
    fdcf:	89 45 fc             	mov    %eax,-0x4(%rbp)
    fdd2:	8b 45 fc             	mov    -0x4(%rbp),%eax
    fdd5:	5d                   	pop    %rbp
    fdd6:	c3                   	ret    
    fdd7:	90                   	nop

000000000000fdd8 <_ZN3MPI6Status7Set_tagEi>:
    fdd8:	f3 0f 1e fa          	endbr64 
    fddc:	55                   	push   %rbp
    fddd:	48 89 e5             	mov    %rsp,%rbp
    fde0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fde4:	89 75 f4             	mov    %esi,-0xc(%rbp)
    fde7:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fdeb:	8b 55 f4             	mov    -0xc(%rbp),%edx
    fdee:	89 50 0c             	mov    %edx,0xc(%rax)
    fdf1:	90                   	nop
    fdf2:	5d                   	pop    %rbp
    fdf3:	c3                   	ret    

000000000000fdf4 <_ZNK3MPI6Status9Get_errorEv>:
    fdf4:	f3 0f 1e fa          	endbr64 
    fdf8:	55                   	push   %rbp
    fdf9:	48 89 e5             	mov    %rsp,%rbp
    fdfc:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fe00:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    fe04:	8b 40 10             	mov    0x10(%rax),%eax
    fe07:	89 45 fc             	mov    %eax,-0x4(%rbp)
    fe0a:	8b 45 fc             	mov    -0x4(%rbp),%eax
    fe0d:	5d                   	pop    %rbp
    fe0e:	c3                   	ret    
    fe0f:	90                   	nop

000000000000fe10 <_ZN3MPI6Status9Set_errorEi>:
    fe10:	f3 0f 1e fa          	endbr64 
    fe14:	55                   	push   %rbp
    fe15:	48 89 e5             	mov    %rsp,%rbp
    fe18:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fe1c:	89 75 f4             	mov    %esi,-0xc(%rbp)
    fe1f:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fe23:	8b 55 f4             	mov    -0xc(%rbp),%edx
    fe26:	89 50 10             	mov    %edx,0x10(%rax)
    fe29:	90                   	nop
    fe2a:	5d                   	pop    %rbp
    fe2b:	c3                   	ret    

000000000000fe2c <_ZN3MPI6Status12Set_elementsERKNS_8DatatypeEi>:
    fe2c:	f3 0f 1e fa          	endbr64 
    fe30:	55                   	push   %rbp
    fe31:	48 89 e5             	mov    %rsp,%rbp
    fe34:	48 83 ec 20          	sub    $0x20,%rsp
    fe38:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fe3c:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    fe40:	89 55 ec             	mov    %edx,-0x14(%rbp)
    fe43:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
    fe47:	48 89 c7             	mov    %rax,%rdi
    fe4a:	e8 51 ba ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
    fe4f:	48 89 c6             	mov    %rax,%rsi
    fe52:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    fe56:	48 8d 48 08          	lea    0x8(%rax),%rcx
    fe5a:	8b 45 ec             	mov    -0x14(%rbp),%eax
    fe5d:	89 c2                	mov    %eax,%edx
    fe5f:	48 89 cf             	mov    %rcx,%rdi
    fe62:	e8 79 a7 ff ff       	call   a5e0 <MPI_Status_set_elements@plt>
    fe67:	90                   	nop
    fe68:	c9                   	leave  
    fe69:	c3                   	ret    

000000000000fe6a <_ZN3MPI6Status13Set_cancelledEb>:
    fe6a:	f3 0f 1e fa          	endbr64 
    fe6e:	55                   	push   %rbp
    fe6f:	48 89 e5             	mov    %rsp,%rbp
    fe72:	48 83 ec 10          	sub    $0x10,%rsp
    fe76:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fe7a:	89 f0                	mov    %esi,%eax
    fe7c:	88 45 f4             	mov    %al,-0xc(%rbp)
    fe7f:	0f b6 45 f4          	movzbl -0xc(%rbp),%eax
    fe83:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    fe87:	48 83 c2 08          	add    $0x8,%rdx
    fe8b:	89 c6                	mov    %eax,%esi
    fe8d:	48 89 d7             	mov    %rdx,%rdi
    fe90:	e8 cb ab ff ff       	call   aa60 <MPI_Status_set_cancelled@plt>
    fe95:	90                   	nop
    fe96:	c9                   	leave  
    fe97:	c3                   	ret    

000000000000fe98 <_ZN3MPI4Info6DeleteEPKc>:
    fe98:	f3 0f 1e fa          	endbr64 
    fe9c:	55                   	push   %rbp
    fe9d:	48 89 e5             	mov    %rsp,%rbp
    fea0:	48 83 ec 10          	sub    $0x10,%rsp
    fea4:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    fea8:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
    feac:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    feb0:	48 8b 40 08          	mov    0x8(%rax),%rax
    feb4:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    feb8:	48 89 d6             	mov    %rdx,%rsi
    febb:	48 89 c7             	mov    %rax,%rdi
    febe:	e8 5d a8 ff ff       	call   a720 <MPI_Info_delete@plt>
    fec3:	90                   	nop
    fec4:	c9                   	leave  
    fec5:	c3                   	ret    

000000000000fec6 <_ZNK3MPI4Info3DupEv>:
    fec6:	f3 0f 1e fa          	endbr64 
    feca:	55                   	push   %rbp
    fecb:	48 89 e5             	mov    %rsp,%rbp
    fece:	48 83 ec 20          	sub    $0x20,%rsp
    fed2:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    fed6:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    feda:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    fee1:	00 00 
    fee3:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    fee7:	31 c0                	xor    %eax,%eax
    fee9:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
    feed:	48 8b 40 08          	mov    0x8(%rax),%rax
    fef1:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
    fef5:	48 89 d6             	mov    %rdx,%rsi
    fef8:	48 89 c7             	mov    %rax,%rdi
    fefb:	e8 30 9f ff ff       	call   9e30 <MPI_Info_dup@plt>
    ff00:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
    ff04:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ff08:	48 89 d6             	mov    %rdx,%rsi
    ff0b:	48 89 c7             	mov    %rax,%rdi
    ff0e:	e8 f9 be ff ff       	call   be0c <_ZN3MPI4InfoC1EP11ompi_info_t>
    ff13:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ff17:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
    ff1e:	00 00 
    ff20:	74 05                	je     ff27 <_ZNK3MPI4Info3DupEv+0x61>
    ff22:	e8 a9 a5 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    ff27:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ff2b:	c9                   	leave  
    ff2c:	c3                   	ret    
    ff2d:	90                   	nop

000000000000ff2e <_ZN3MPI4Info4FreeEv>:
    ff2e:	f3 0f 1e fa          	endbr64 
    ff32:	55                   	push   %rbp
    ff33:	48 89 e5             	mov    %rsp,%rbp
    ff36:	48 83 ec 10          	sub    $0x10,%rsp
    ff3a:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
    ff3e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
    ff42:	48 83 c0 08          	add    $0x8,%rax
    ff46:	48 89 c7             	mov    %rax,%rdi
    ff49:	e8 32 ab ff ff       	call   aa80 <MPI_Info_free@plt>
    ff4e:	90                   	nop
    ff4f:	c9                   	leave  
    ff50:	c3                   	ret    
    ff51:	90                   	nop

000000000000ff52 <_ZNK3MPI4Info3GetEPKciPc>:
    ff52:	f3 0f 1e fa          	endbr64 
    ff56:	55                   	push   %rbp
    ff57:	48 89 e5             	mov    %rsp,%rbp
    ff5a:	48 83 ec 30          	sub    $0x30,%rsp
    ff5e:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    ff62:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
    ff66:	89 55 dc             	mov    %edx,-0x24(%rbp)
    ff69:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
    ff6d:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ff74:	00 00 
    ff76:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    ff7a:	31 c0                	xor    %eax,%eax
    ff7c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ff80:	48 8b 40 08          	mov    0x8(%rax),%rax
    ff84:	48 8d 7d f4          	lea    -0xc(%rbp),%rdi
    ff88:	48 8b 4d d0          	mov    -0x30(%rbp),%rcx
    ff8c:	8b 55 dc             	mov    -0x24(%rbp),%edx
    ff8f:	48 8b 75 e0          	mov    -0x20(%rbp),%rsi
    ff93:	49 89 f8             	mov    %rdi,%r8
    ff96:	48 89 c7             	mov    %rax,%rdi
    ff99:	e8 12 a1 ff ff       	call   a0b0 <MPI_Info_get@plt>
    ff9e:	8b 45 f4             	mov    -0xc(%rbp),%eax
    ffa1:	85 c0                	test   %eax,%eax
    ffa3:	0f 95 c0             	setne  %al
    ffa6:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    ffaa:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
    ffb1:	00 00 
    ffb3:	74 05                	je     ffba <_ZNK3MPI4Info3GetEPKciPc+0x68>
    ffb5:	e8 16 a5 ff ff       	call   a4d0 <__stack_chk_fail@plt>
    ffba:	c9                   	leave  
    ffbb:	c3                   	ret    

000000000000ffbc <_ZNK3MPI4Info9Get_nkeysEv>:
    ffbc:	f3 0f 1e fa          	endbr64 
    ffc0:	55                   	push   %rbp
    ffc1:	48 89 e5             	mov    %rsp,%rbp
    ffc4:	48 83 ec 20          	sub    $0x20,%rsp
    ffc8:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
    ffcc:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
    ffd3:	00 00 
    ffd5:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
    ffd9:	31 c0                	xor    %eax,%eax
    ffdb:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
    ffdf:	48 8b 40 08          	mov    0x8(%rax),%rax
    ffe3:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
    ffe7:	48 89 d6             	mov    %rdx,%rsi
    ffea:	48 89 c7             	mov    %rax,%rdi
    ffed:	e8 ae a4 ff ff       	call   a4a0 <MPI_Info_get_nkeys@plt>
    fff2:	8b 45 f4             	mov    -0xc(%rbp),%eax
    fff5:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
    fff9:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
   10000:	00 00 
   10002:	74 05                	je     10009 <_ZNK3MPI4Info9Get_nkeysEv+0x4d>
   10004:	e8 c7 a4 ff ff       	call   a4d0 <__stack_chk_fail@plt>
   10009:	c9                   	leave  
   1000a:	c3                   	ret    
   1000b:	90                   	nop

000000000001000c <_ZNK3MPI4Info10Get_nthkeyEiPc>:
   1000c:	f3 0f 1e fa          	endbr64 
   10010:	55                   	push   %rbp
   10011:	48 89 e5             	mov    %rsp,%rbp
   10014:	48 83 ec 20          	sub    $0x20,%rsp
   10018:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   1001c:	89 75 f4             	mov    %esi,-0xc(%rbp)
   1001f:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
   10023:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10027:	48 8b 40 08          	mov    0x8(%rax),%rax
   1002b:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   1002f:	8b 4d f4             	mov    -0xc(%rbp),%ecx
   10032:	89 ce                	mov    %ecx,%esi
   10034:	48 89 c7             	mov    %rax,%rdi
   10037:	e8 24 a7 ff ff       	call   a760 <MPI_Info_get_nthkey@plt>
   1003c:	90                   	nop
   1003d:	c9                   	leave  
   1003e:	c3                   	ret    
   1003f:	90                   	nop

0000000000010040 <_ZNK3MPI4Info12Get_valuelenEPKcRi>:
   10040:	f3 0f 1e fa          	endbr64 
   10044:	55                   	push   %rbp
   10045:	48 89 e5             	mov    %rsp,%rbp
   10048:	48 83 ec 30          	sub    $0x30,%rsp
   1004c:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   10050:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   10054:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
   10058:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   1005f:	00 00 
   10061:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   10065:	31 c0                	xor    %eax,%eax
   10067:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1006b:	48 8b 40 08          	mov    0x8(%rax),%rax
   1006f:	48 8d 4d f4          	lea    -0xc(%rbp),%rcx
   10073:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
   10077:	48 8b 75 e0          	mov    -0x20(%rbp),%rsi
   1007b:	48 89 c7             	mov    %rax,%rdi
   1007e:	e8 3d a8 ff ff       	call   a8c0 <MPI_Info_get_valuelen@plt>
   10083:	8b 45 f4             	mov    -0xc(%rbp),%eax
   10086:	85 c0                	test   %eax,%eax
   10088:	0f 95 c0             	setne  %al
   1008b:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   1008f:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
   10096:	00 00 
   10098:	74 05                	je     1009f <_ZNK3MPI4Info12Get_valuelenEPKcRi+0x5f>
   1009a:	e8 31 a4 ff ff       	call   a4d0 <__stack_chk_fail@plt>
   1009f:	c9                   	leave  
   100a0:	c3                   	ret    
   100a1:	90                   	nop

00000000000100a2 <_ZN3MPI4Info3SetEPKcS2_>:
   100a2:	f3 0f 1e fa          	endbr64 
   100a6:	55                   	push   %rbp
   100a7:	48 89 e5             	mov    %rsp,%rbp
   100aa:	48 83 ec 20          	sub    $0x20,%rsp
   100ae:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   100b2:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   100b6:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
   100ba:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   100be:	48 8b 40 08          	mov    0x8(%rax),%rax
   100c2:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   100c6:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
   100ca:	48 89 ce             	mov    %rcx,%rsi
   100cd:	48 89 c7             	mov    %rax,%rdi
   100d0:	e8 db 9e ff ff       	call   9fb0 <MPI_Info_set@plt>
   100d5:	90                   	nop
   100d6:	c9                   	leave  
   100d7:	c3                   	ret    

00000000000100d8 <_ZNK3MPI3Win14Get_errhandlerEv>:
   100d8:	f3 0f 1e fa          	endbr64 
   100dc:	55                   	push   %rbp
   100dd:	48 89 e5             	mov    %rsp,%rbp
   100e0:	48 83 ec 20          	sub    $0x20,%rsp
   100e4:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   100e8:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   100ec:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   100f3:	00 00 
   100f5:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   100f9:	31 c0                	xor    %eax,%eax
   100fb:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   100ff:	48 8b 40 08          	mov    0x8(%rax),%rax
   10103:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
   10107:	48 89 d6             	mov    %rdx,%rsi
   1010a:	48 89 c7             	mov    %rax,%rdi
   1010d:	e8 0e a9 ff ff       	call   aa20 <MPI_Win_get_errhandler@plt>
   10112:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   10116:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1011a:	48 89 d6             	mov    %rdx,%rsi
   1011d:	48 89 c7             	mov    %rax,%rdi
   10120:	e8 81 bb ff ff       	call   bca6 <_ZN3MPI10ErrhandlerC1EP17ompi_errhandler_t>
   10125:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10129:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   10130:	00 00 
   10132:	74 05                	je     10139 <_ZNK3MPI3Win14Get_errhandlerEv+0x61>
   10134:	e8 97 a3 ff ff       	call   a4d0 <__stack_chk_fail@plt>
   10139:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1013d:	c9                   	leave  
   1013e:	c3                   	ret    
   1013f:	90                   	nop

0000000000010140 <_ZNK3MPI3Win14Set_errhandlerERKNS_10ErrhandlerE>:
   10140:	f3 0f 1e fa          	endbr64 
   10144:	55                   	push   %rbp
   10145:	48 89 e5             	mov    %rsp,%rbp
   10148:	48 83 ec 10          	sub    $0x10,%rsp
   1014c:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10150:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   10154:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   10158:	48 89 c7             	mov    %rax,%rdi
   1015b:	e8 74 bb ff ff       	call   bcd4 <_ZNK3MPI10ErrhandlercvP17ompi_errhandler_tEv>
   10160:	48 89 c2             	mov    %rax,%rdx
   10163:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10167:	48 8b 40 08          	mov    0x8(%rax),%rax
   1016b:	48 89 d6             	mov    %rdx,%rsi
   1016e:	48 89 c7             	mov    %rax,%rdi
   10171:	e8 0a a6 ff ff       	call   a780 <MPI_Win_set_errhandler@plt>
   10176:	90                   	nop
   10177:	c9                   	leave  
   10178:	c3                   	ret    
   10179:	90                   	nop

000000000001017a <_ZNK3MPI3Win10AccumulateEPKviRKNS_8DatatypeEiliS5_RKNS_2OpE>:
   1017a:	f3 0f 1e fa          	endbr64 
   1017e:	55                   	push   %rbp
   1017f:	48 89 e5             	mov    %rsp,%rbp
   10182:	41 55                	push   %r13
   10184:	41 54                	push   %r12
   10186:	53                   	push   %rbx
   10187:	48 83 ec 38          	sub    $0x38,%rsp
   1018b:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   1018f:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
   10193:	89 55 cc             	mov    %edx,-0x34(%rbp)
   10196:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
   1019a:	44 89 45 c8          	mov    %r8d,-0x38(%rbp)
   1019e:	4c 89 4d b8          	mov    %r9,-0x48(%rbp)
   101a2:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   101a6:	4c 8b 68 08          	mov    0x8(%rax),%r13
   101aa:	48 8b 45 20          	mov    0x20(%rbp),%rax
   101ae:	48 89 c7             	mov    %rax,%rdi
   101b1:	e8 60 f9 ff ff       	call   fb16 <_ZNK3MPI2OpcvP9ompi_op_tEv>
   101b6:	49 89 c4             	mov    %rax,%r12
   101b9:	48 8b 45 18          	mov    0x18(%rbp),%rax
   101bd:	48 89 c7             	mov    %rax,%rdi
   101c0:	e8 db b6 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
   101c5:	48 89 c3             	mov    %rax,%rbx
   101c8:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   101cc:	48 89 c7             	mov    %rax,%rdi
   101cf:	e8 cc b6 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
   101d4:	48 89 c7             	mov    %rax,%rdi
   101d7:	48 8b 4d b8          	mov    -0x48(%rbp),%rcx
   101db:	8b 55 c8             	mov    -0x38(%rbp),%edx
   101de:	8b 75 cc             	mov    -0x34(%rbp),%esi
   101e1:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   101e5:	48 83 ec 08          	sub    $0x8,%rsp
   101e9:	41 55                	push   %r13
   101eb:	41 54                	push   %r12
   101ed:	53                   	push   %rbx
   101ee:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
   101f2:	49 89 c8             	mov    %rcx,%r8
   101f5:	89 d1                	mov    %edx,%ecx
   101f7:	48 89 fa             	mov    %rdi,%rdx
   101fa:	48 89 c7             	mov    %rax,%rdi
   101fd:	e8 2e a8 ff ff       	call   aa30 <MPI_Accumulate@plt>
   10202:	48 83 c4 20          	add    $0x20,%rsp
   10206:	90                   	nop
   10207:	48 8d 65 e8          	lea    -0x18(%rbp),%rsp
   1020b:	5b                   	pop    %rbx
   1020c:	41 5c                	pop    %r12
   1020e:	41 5d                	pop    %r13
   10210:	5d                   	pop    %rbp
   10211:	c3                   	ret    

0000000000010212 <_ZNK3MPI3Win8CompleteEv>:
   10212:	f3 0f 1e fa          	endbr64 
   10216:	55                   	push   %rbp
   10217:	48 89 e5             	mov    %rsp,%rbp
   1021a:	48 83 ec 10          	sub    $0x10,%rsp
   1021e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10222:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10226:	48 8b 40 08          	mov    0x8(%rax),%rax
   1022a:	48 89 c7             	mov    %rax,%rdi
   1022d:	e8 0e a7 ff ff       	call   a940 <MPI_Win_complete@plt>
   10232:	90                   	nop
   10233:	c9                   	leave  
   10234:	c3                   	ret    
   10235:	90                   	nop

0000000000010236 <_ZNK3MPI3Win5FenceEi>:
   10236:	f3 0f 1e fa          	endbr64 
   1023a:	55                   	push   %rbp
   1023b:	48 89 e5             	mov    %rsp,%rbp
   1023e:	48 83 ec 10          	sub    $0x10,%rsp
   10242:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10246:	89 75 f4             	mov    %esi,-0xc(%rbp)
   10249:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1024d:	48 8b 50 08          	mov    0x8(%rax),%rdx
   10251:	8b 45 f4             	mov    -0xc(%rbp),%eax
   10254:	48 89 d6             	mov    %rdx,%rsi
   10257:	89 c7                	mov    %eax,%edi
   10259:	e8 d2 a6 ff ff       	call   a930 <MPI_Win_fence@plt>
   1025e:	90                   	nop
   1025f:	c9                   	leave  
   10260:	c3                   	ret    
   10261:	90                   	nop

0000000000010262 <_ZNK3MPI3Win3GetEPKviRKNS_8DatatypeEiliS5_>:
   10262:	f3 0f 1e fa          	endbr64 
   10266:	55                   	push   %rbp
   10267:	48 89 e5             	mov    %rsp,%rbp
   1026a:	41 54                	push   %r12
   1026c:	53                   	push   %rbx
   1026d:	48 83 ec 30          	sub    $0x30,%rsp
   10271:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   10275:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   10279:	89 55 dc             	mov    %edx,-0x24(%rbp)
   1027c:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
   10280:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
   10284:	4c 89 4d c8          	mov    %r9,-0x38(%rbp)
   10288:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1028c:	4c 8b 60 08          	mov    0x8(%rax),%r12
   10290:	48 8b 45 18          	mov    0x18(%rbp),%rax
   10294:	48 89 c7             	mov    %rax,%rdi
   10297:	e8 04 b6 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
   1029c:	48 89 c3             	mov    %rax,%rbx
   1029f:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   102a3:	48 89 c7             	mov    %rax,%rdi
   102a6:	e8 f5 b5 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
   102ab:	48 89 c7             	mov    %rax,%rdi
   102ae:	48 8b 4d c8          	mov    -0x38(%rbp),%rcx
   102b2:	8b 55 d8             	mov    -0x28(%rbp),%edx
   102b5:	8b 75 dc             	mov    -0x24(%rbp),%esi
   102b8:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   102bc:	41 54                	push   %r12
   102be:	53                   	push   %rbx
   102bf:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
   102c3:	49 89 c8             	mov    %rcx,%r8
   102c6:	89 d1                	mov    %edx,%ecx
   102c8:	48 89 fa             	mov    %rdi,%rdx
   102cb:	48 89 c7             	mov    %rax,%rdi
   102ce:	e8 1d 9e ff ff       	call   a0f0 <MPI_Get@plt>
   102d3:	48 83 c4 10          	add    $0x10,%rsp
   102d7:	90                   	nop
   102d8:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
   102dc:	5b                   	pop    %rbx
   102dd:	41 5c                	pop    %r12
   102df:	5d                   	pop    %rbp
   102e0:	c3                   	ret    
   102e1:	90                   	nop

00000000000102e2 <_ZNK3MPI3Win9Get_groupEv>:
   102e2:	f3 0f 1e fa          	endbr64 
   102e6:	55                   	push   %rbp
   102e7:	48 89 e5             	mov    %rsp,%rbp
   102ea:	48 83 ec 20          	sub    $0x20,%rsp
   102ee:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   102f2:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   102f6:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   102fd:	00 00 
   102ff:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   10303:	31 c0                	xor    %eax,%eax
   10305:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   10309:	48 8b 40 08          	mov    0x8(%rax),%rax
   1030d:	48 8d 55 f0          	lea    -0x10(%rbp),%rdx
   10311:	48 89 d6             	mov    %rdx,%rsi
   10314:	48 89 c7             	mov    %rax,%rdi
   10317:	e8 44 9e ff ff       	call   a160 <MPI_Win_get_group@plt>
   1031c:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   10320:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10324:	48 89 d6             	mov    %rdx,%rsi
   10327:	48 89 c7             	mov    %rax,%rdi
   1032a:	e8 7f b7 ff ff       	call   baae <_ZN3MPI5GroupC1EP12ompi_group_t>
   1032f:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10333:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   1033a:	00 00 
   1033c:	74 05                	je     10343 <_ZNK3MPI3Win9Get_groupEv+0x61>
   1033e:	e8 8d a1 ff ff       	call   a4d0 <__stack_chk_fail@plt>
   10343:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10347:	c9                   	leave  
   10348:	c3                   	ret    
   10349:	90                   	nop

000000000001034a <_ZNK3MPI3Win4LockEiii>:
   1034a:	f3 0f 1e fa          	endbr64 
   1034e:	55                   	push   %rbp
   1034f:	48 89 e5             	mov    %rsp,%rbp
   10352:	48 83 ec 20          	sub    $0x20,%rsp
   10356:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   1035a:	89 75 f4             	mov    %esi,-0xc(%rbp)
   1035d:	89 55 f0             	mov    %edx,-0x10(%rbp)
   10360:	89 4d ec             	mov    %ecx,-0x14(%rbp)
   10363:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10367:	48 8b 48 08          	mov    0x8(%rax),%rcx
   1036b:	8b 55 ec             	mov    -0x14(%rbp),%edx
   1036e:	8b 75 f0             	mov    -0x10(%rbp),%esi
   10371:	8b 45 f4             	mov    -0xc(%rbp),%eax
   10374:	89 c7                	mov    %eax,%edi
   10376:	e8 c5 a1 ff ff       	call   a540 <MPI_Win_lock@plt>
   1037b:	90                   	nop
   1037c:	c9                   	leave  
   1037d:	c3                   	ret    

000000000001037e <_ZNK3MPI3Win4PostERKNS_5GroupEi>:
   1037e:	f3 0f 1e fa          	endbr64 
   10382:	55                   	push   %rbp
   10383:	48 89 e5             	mov    %rsp,%rbp
   10386:	53                   	push   %rbx
   10387:	48 83 ec 28          	sub    $0x28,%rsp
   1038b:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   1038f:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   10393:	89 55 dc             	mov    %edx,-0x24(%rbp)
   10396:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1039a:	48 8b 58 08          	mov    0x8(%rax),%rbx
   1039e:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   103a2:	48 89 c7             	mov    %rax,%rdi
   103a5:	e8 80 b7 ff ff       	call   bb2a <_ZNK3MPI5GroupcvP12ompi_group_tEv>
   103aa:	48 89 c1             	mov    %rax,%rcx
   103ad:	8b 45 dc             	mov    -0x24(%rbp),%eax
   103b0:	48 89 da             	mov    %rbx,%rdx
   103b3:	89 c6                	mov    %eax,%esi
   103b5:	48 89 cf             	mov    %rcx,%rdi
   103b8:	e8 23 9c ff ff       	call   9fe0 <MPI_Win_post@plt>
   103bd:	90                   	nop
   103be:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   103c2:	c9                   	leave  
   103c3:	c3                   	ret    

00000000000103c4 <_ZNK3MPI3Win3PutEPKviRKNS_8DatatypeEiliS5_>:
   103c4:	f3 0f 1e fa          	endbr64 
   103c8:	55                   	push   %rbp
   103c9:	48 89 e5             	mov    %rsp,%rbp
   103cc:	41 54                	push   %r12
   103ce:	53                   	push   %rbx
   103cf:	48 83 ec 30          	sub    $0x30,%rsp
   103d3:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   103d7:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   103db:	89 55 dc             	mov    %edx,-0x24(%rbp)
   103de:	48 89 4d d0          	mov    %rcx,-0x30(%rbp)
   103e2:	44 89 45 d8          	mov    %r8d,-0x28(%rbp)
   103e6:	4c 89 4d c8          	mov    %r9,-0x38(%rbp)
   103ea:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   103ee:	4c 8b 60 08          	mov    0x8(%rax),%r12
   103f2:	48 8b 45 18          	mov    0x18(%rbp),%rax
   103f6:	48 89 c7             	mov    %rax,%rdi
   103f9:	e8 a2 b4 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
   103fe:	48 89 c3             	mov    %rax,%rbx
   10401:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   10405:	48 89 c7             	mov    %rax,%rdi
   10408:	e8 93 b4 ff ff       	call   b8a0 <_ZNK3MPI8DatatypecvP15ompi_datatype_tEv>
   1040d:	48 89 c7             	mov    %rax,%rdi
   10410:	48 8b 4d c8          	mov    -0x38(%rbp),%rcx
   10414:	8b 55 d8             	mov    -0x28(%rbp),%edx
   10417:	8b 75 dc             	mov    -0x24(%rbp),%esi
   1041a:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   1041e:	41 54                	push   %r12
   10420:	53                   	push   %rbx
   10421:	44 8b 4d 10          	mov    0x10(%rbp),%r9d
   10425:	49 89 c8             	mov    %rcx,%r8
   10428:	89 d1                	mov    %edx,%ecx
   1042a:	48 89 fa             	mov    %rdi,%rdx
   1042d:	48 89 c7             	mov    %rax,%rdi
   10430:	e8 7b a1 ff ff       	call   a5b0 <MPI_Put@plt>
   10435:	48 83 c4 10          	add    $0x10,%rsp
   10439:	90                   	nop
   1043a:	48 8d 65 f0          	lea    -0x10(%rbp),%rsp
   1043e:	5b                   	pop    %rbx
   1043f:	41 5c                	pop    %r12
   10441:	5d                   	pop    %rbp
   10442:	c3                   	ret    
   10443:	90                   	nop

0000000000010444 <_ZNK3MPI3Win5StartERKNS_5GroupEi>:
   10444:	f3 0f 1e fa          	endbr64 
   10448:	55                   	push   %rbp
   10449:	48 89 e5             	mov    %rsp,%rbp
   1044c:	53                   	push   %rbx
   1044d:	48 83 ec 28          	sub    $0x28,%rsp
   10451:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   10455:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   10459:	89 55 dc             	mov    %edx,-0x24(%rbp)
   1045c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10460:	48 8b 58 08          	mov    0x8(%rax),%rbx
   10464:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   10468:	48 89 c7             	mov    %rax,%rdi
   1046b:	e8 ba b6 ff ff       	call   bb2a <_ZNK3MPI5GroupcvP12ompi_group_tEv>
   10470:	48 89 c1             	mov    %rax,%rcx
   10473:	8b 45 dc             	mov    -0x24(%rbp),%eax
   10476:	48 89 da             	mov    %rbx,%rdx
   10479:	89 c6                	mov    %eax,%esi
   1047b:	48 89 cf             	mov    %rcx,%rdi
   1047e:	e8 4d 9d ff ff       	call   a1d0 <MPI_Win_start@plt>
   10483:	90                   	nop
   10484:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   10488:	c9                   	leave  
   10489:	c3                   	ret    

000000000001048a <_ZNK3MPI3Win4TestEv>:
   1048a:	f3 0f 1e fa          	endbr64 
   1048e:	55                   	push   %rbp
   1048f:	48 89 e5             	mov    %rsp,%rbp
   10492:	48 83 ec 20          	sub    $0x20,%rsp
   10496:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   1049a:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   104a1:	00 00 
   104a3:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   104a7:	31 c0                	xor    %eax,%eax
   104a9:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   104ad:	48 8b 40 08          	mov    0x8(%rax),%rax
   104b1:	48 8d 55 f4          	lea    -0xc(%rbp),%rdx
   104b5:	48 89 d6             	mov    %rdx,%rsi
   104b8:	48 89 c7             	mov    %rax,%rdi
   104bb:	e8 10 a4 ff ff       	call   a8d0 <MPI_Win_test@plt>
   104c0:	8b 45 f4             	mov    -0xc(%rbp),%eax
   104c3:	85 c0                	test   %eax,%eax
   104c5:	0f 95 c0             	setne  %al
   104c8:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   104cc:	64 48 2b 14 25 28 00 	sub    %fs:0x28,%rdx
   104d3:	00 00 
   104d5:	74 05                	je     104dc <_ZNK3MPI3Win4TestEv+0x52>
   104d7:	e8 f4 9f ff ff       	call   a4d0 <__stack_chk_fail@plt>
   104dc:	c9                   	leave  
   104dd:	c3                   	ret    

00000000000104de <_ZNK3MPI3Win6UnlockEi>:
   104de:	f3 0f 1e fa          	endbr64 
   104e2:	55                   	push   %rbp
   104e3:	48 89 e5             	mov    %rsp,%rbp
   104e6:	48 83 ec 10          	sub    $0x10,%rsp
   104ea:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   104ee:	89 75 f4             	mov    %esi,-0xc(%rbp)
   104f1:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   104f5:	48 8b 50 08          	mov    0x8(%rax),%rdx
   104f9:	8b 45 f4             	mov    -0xc(%rbp),%eax
   104fc:	48 89 d6             	mov    %rdx,%rsi
   104ff:	89 c7                	mov    %eax,%edi
   10501:	e8 ca 98 ff ff       	call   9dd0 <MPI_Win_unlock@plt>
   10506:	90                   	nop
   10507:	c9                   	leave  
   10508:	c3                   	ret    
   10509:	90                   	nop

000000000001050a <_ZNK3MPI3Win4WaitEv>:
   1050a:	f3 0f 1e fa          	endbr64 
   1050e:	55                   	push   %rbp
   1050f:	48 89 e5             	mov    %rsp,%rbp
   10512:	48 83 ec 10          	sub    $0x10,%rsp
   10516:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   1051a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1051e:	48 8b 40 08          	mov    0x8(%rax),%rax
   10522:	48 89 c7             	mov    %rax,%rdi
   10525:	e8 86 9d ff ff       	call   a2b0 <MPI_Win_wait@plt>
   1052a:	90                   	nop
   1052b:	c9                   	leave  
   1052c:	c3                   	ret    
   1052d:	90                   	nop

000000000001052e <_ZNK3MPI3Win15Call_errhandlerEi>:
   1052e:	f3 0f 1e fa          	endbr64 
   10532:	55                   	push   %rbp
   10533:	48 89 e5             	mov    %rsp,%rbp
   10536:	48 83 ec 10          	sub    $0x10,%rsp
   1053a:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   1053e:	89 75 f4             	mov    %esi,-0xc(%rbp)
   10541:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10545:	48 8b 40 08          	mov    0x8(%rax),%rax
   10549:	8b 55 f4             	mov    -0xc(%rbp),%edx
   1054c:	89 d6                	mov    %edx,%esi
   1054e:	48 89 c7             	mov    %rax,%rdi
   10551:	e8 aa a3 ff ff       	call   a900 <MPI_Win_call_errhandler@plt>
   10556:	90                   	nop
   10557:	c9                   	leave  
   10558:	c3                   	ret    
   10559:	90                   	nop

000000000001055a <_ZN3MPI3Win11Delete_attrEi>:
   1055a:	f3 0f 1e fa          	endbr64 
   1055e:	55                   	push   %rbp
   1055f:	48 89 e5             	mov    %rsp,%rbp
   10562:	48 83 ec 10          	sub    $0x10,%rsp
   10566:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   1056a:	89 75 f4             	mov    %esi,-0xc(%rbp)
   1056d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10571:	48 8b 40 08          	mov    0x8(%rax),%rax
   10575:	8b 55 f4             	mov    -0xc(%rbp),%edx
   10578:	89 d6                	mov    %edx,%esi
   1057a:	48 89 c7             	mov    %rax,%rdi
   1057d:	e8 4e a4 ff ff       	call   a9d0 <MPI_Win_delete_attr@plt>
   10582:	90                   	nop
   10583:	c9                   	leave  
   10584:	c3                   	ret    
   10585:	90                   	nop

0000000000010586 <_ZNK3MPI3Win8Get_nameEPcRi>:
   10586:	f3 0f 1e fa          	endbr64 
   1058a:	55                   	push   %rbp
   1058b:	48 89 e5             	mov    %rsp,%rbp
   1058e:	48 83 ec 20          	sub    $0x20,%rsp
   10592:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10596:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   1059a:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
   1059e:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   105a2:	48 8b 40 08          	mov    0x8(%rax),%rax
   105a6:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   105aa:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
   105ae:	48 89 ce             	mov    %rcx,%rsi
   105b1:	48 89 c7             	mov    %rax,%rdi
   105b4:	e8 37 99 ff ff       	call   9ef0 <MPI_Win_get_name@plt>
   105b9:	90                   	nop
   105ba:	c9                   	leave  
   105bb:	c3                   	ret    

00000000000105bc <_ZN3MPI3Win8Set_attrEiPKv>:
   105bc:	f3 0f 1e fa          	endbr64 
   105c0:	55                   	push   %rbp
   105c1:	48 89 e5             	mov    %rsp,%rbp
   105c4:	48 83 ec 20          	sub    $0x20,%rsp
   105c8:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   105cc:	89 75 f4             	mov    %esi,-0xc(%rbp)
   105cf:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
   105d3:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   105d7:	48 8b 40 08          	mov    0x8(%rax),%rax
   105db:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   105df:	8b 4d f4             	mov    -0xc(%rbp),%ecx
   105e2:	89 ce                	mov    %ecx,%esi
   105e4:	48 89 c7             	mov    %rax,%rdi
   105e7:	e8 d4 a3 ff ff       	call   a9c0 <MPI_Win_set_attr@plt>
   105ec:	90                   	nop
   105ed:	c9                   	leave  
   105ee:	c3                   	ret    
   105ef:	90                   	nop

00000000000105f0 <_ZN3MPI3Win8Set_nameEPKc>:
   105f0:	f3 0f 1e fa          	endbr64 
   105f4:	55                   	push   %rbp
   105f5:	48 89 e5             	mov    %rsp,%rbp
   105f8:	48 83 ec 10          	sub    $0x10,%rsp
   105fc:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10600:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   10604:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10608:	48 8b 40 08          	mov    0x8(%rax),%rax
   1060c:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   10610:	48 89 d6             	mov    %rdx,%rsi
   10613:	48 89 c7             	mov    %rax,%rdi
   10616:	e8 45 9d ff ff       	call   a360 <MPI_Win_set_name@plt>
   1061b:	90                   	nop
   1061c:	c9                   	leave  
   1061d:	c3                   	ret    

000000000001061e <_ZNSt8__detail18__to_chars_10_implIjEEvPcjT_>:
   1061e:	f3 0f 1e fa          	endbr64 
   10622:	55                   	push   %rbp
   10623:	48 89 e5             	mov    %rsp,%rbp
   10626:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   1062a:	89 75 e4             	mov    %esi,-0x1c(%rbp)
   1062d:	89 55 e0             	mov    %edx,-0x20(%rbp)
   10630:	8b 45 e4             	mov    -0x1c(%rbp),%eax
   10633:	83 e8 01             	sub    $0x1,%eax
   10636:	89 45 f4             	mov    %eax,-0xc(%rbp)
   10639:	eb 77                	jmp    106b2 <_ZNSt8__detail18__to_chars_10_implIjEEvPcjT_+0x94>
   1063b:	8b 55 e0             	mov    -0x20(%rbp),%edx
   1063e:	89 d0                	mov    %edx,%eax
   10640:	48 69 c0 1f 85 eb 51 	imul   $0x51eb851f,%rax,%rax
   10647:	48 c1 e8 20          	shr    $0x20,%rax
   1064b:	c1 e8 05             	shr    $0x5,%eax
   1064e:	6b c8 64             	imul   $0x64,%eax,%ecx
   10651:	89 d0                	mov    %edx,%eax
   10653:	29 c8                	sub    %ecx,%eax
   10655:	01 c0                	add    %eax,%eax
   10657:	89 45 fc             	mov    %eax,-0x4(%rbp)
   1065a:	8b 45 e0             	mov    -0x20(%rbp),%eax
   1065d:	89 c0                	mov    %eax,%eax
   1065f:	48 69 c0 1f 85 eb 51 	imul   $0x51eb851f,%rax,%rax
   10666:	48 c1 e8 20          	shr    $0x20,%rax
   1066a:	c1 e8 05             	shr    $0x5,%eax
   1066d:	89 45 e0             	mov    %eax,-0x20(%rbp)
   10670:	8b 45 fc             	mov    -0x4(%rbp),%eax
   10673:	8d 48 01             	lea    0x1(%rax),%ecx
   10676:	8b 55 f4             	mov    -0xc(%rbp),%edx
   10679:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1067d:	48 01 c2             	add    %rax,%rdx
   10680:	89 c9                	mov    %ecx,%ecx
   10682:	48 8d 05 d7 5e 00 00 	lea    0x5ed7(%rip),%rax        # 16560 <_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits>
   10689:	0f b6 04 01          	movzbl (%rcx,%rax,1),%eax
   1068d:	88 02                	mov    %al,(%rdx)
   1068f:	8b 45 f4             	mov    -0xc(%rbp),%eax
   10692:	83 e8 01             	sub    $0x1,%eax
   10695:	89 c2                	mov    %eax,%edx
   10697:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1069b:	48 01 c2             	add    %rax,%rdx
   1069e:	8b 45 fc             	mov    -0x4(%rbp),%eax
   106a1:	48 8d 0d b8 5e 00 00 	lea    0x5eb8(%rip),%rcx        # 16560 <_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits>
   106a8:	0f b6 04 08          	movzbl (%rax,%rcx,1),%eax
   106ac:	88 02                	mov    %al,(%rdx)
   106ae:	83 6d f4 02          	subl   $0x2,-0xc(%rbp)
   106b2:	83 7d e0 63          	cmpl   $0x63,-0x20(%rbp)
   106b6:	77 83                	ja     1063b <_ZNSt8__detail18__to_chars_10_implIjEEvPcjT_+0x1d>
   106b8:	83 7d e0 09          	cmpl   $0x9,-0x20(%rbp)
   106bc:	76 3b                	jbe    106f9 <_ZNSt8__detail18__to_chars_10_implIjEEvPcjT_+0xdb>
   106be:	8b 45 e0             	mov    -0x20(%rbp),%eax
   106c1:	01 c0                	add    %eax,%eax
   106c3:	89 45 f8             	mov    %eax,-0x8(%rbp)
   106c6:	8b 45 f8             	mov    -0x8(%rbp),%eax
   106c9:	8d 48 01             	lea    0x1(%rax),%ecx
   106cc:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   106d0:	48 8d 50 01          	lea    0x1(%rax),%rdx
   106d4:	89 c9                	mov    %ecx,%ecx
   106d6:	48 8d 05 83 5e 00 00 	lea    0x5e83(%rip),%rax        # 16560 <_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits>
   106dd:	0f b6 04 01          	movzbl (%rcx,%rax,1),%eax
   106e1:	88 02                	mov    %al,(%rdx)
   106e3:	8b 45 f8             	mov    -0x8(%rbp),%eax
   106e6:	48 8d 15 73 5e 00 00 	lea    0x5e73(%rip),%rdx        # 16560 <_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits>
   106ed:	0f b6 14 10          	movzbl (%rax,%rdx,1),%edx
   106f1:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   106f5:	88 10                	mov    %dl,(%rax)
   106f7:	eb 0e                	jmp    10707 <_ZNSt8__detail18__to_chars_10_implIjEEvPcjT_+0xe9>
   106f9:	8b 45 e0             	mov    -0x20(%rbp),%eax
   106fc:	83 c0 30             	add    $0x30,%eax
   106ff:	89 c2                	mov    %eax,%edx
   10701:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10705:	88 10                	mov    %dl,(%rax)
   10707:	90                   	nop
   10708:	5d                   	pop    %rbp
   10709:	c3                   	ret    

000000000001070a <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_>:
   1070a:	f3 0f 1e fa          	endbr64 
   1070e:	55                   	push   %rbp
   1070f:	48 89 e5             	mov    %rsp,%rbp
   10712:	53                   	push   %rbx
   10713:	48 83 ec 48          	sub    $0x48,%rsp
   10717:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
   1071b:	48 89 75 c0          	mov    %rsi,-0x40(%rbp)
   1071f:	48 89 55 b8          	mov    %rdx,-0x48(%rbp)
   10723:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   1072a:	00 00 
   1072c:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   10730:	31 c0                	xor    %eax,%eax
   10732:	c6 45 df 00          	movb   $0x0,-0x21(%rbp)
   10736:	48 8d 45 de          	lea    -0x22(%rbp),%rax
   1073a:	48 89 c7             	mov    %rax,%rdi
   1073d:	e8 a8 b0 ff ff       	call   b7ea <_ZNKSt17integral_constantIbLb1EEcvbEv>
   10742:	84 c0                	test   %al,%al
   10744:	74 06                	je     1074c <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x42>
   10746:	c6 45 df 01          	movb   $0x1,-0x21(%rbp)
   1074a:	eb 5b                	jmp    107a7 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x9d>
   1074c:	48 8d 45 de          	lea    -0x22(%rbp),%rax
   10750:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
   10754:	48 89 d6             	mov    %rdx,%rsi
   10757:	48 89 c7             	mov    %rax,%rdi
   1075a:	e8 31 98 ff ff       	call   9f90 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE13get_allocatorEv@plt>
   1075f:	48 8d 45 dd          	lea    -0x23(%rbp),%rax
   10763:	48 8b 55 c0          	mov    -0x40(%rbp),%rdx
   10767:	48 89 d6             	mov    %rdx,%rsi
   1076a:	48 89 c7             	mov    %rax,%rdi
   1076d:	e8 1e 98 ff ff       	call   9f90 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE13get_allocatorEv@plt>
   10772:	48 8d 55 de          	lea    -0x22(%rbp),%rdx
   10776:	48 8d 45 dd          	lea    -0x23(%rbp),%rax
   1077a:	48 89 d6             	mov    %rdx,%rsi
   1077d:	48 89 c7             	mov    %rax,%rdi
   10780:	e8 4b 01 00 00       	call   108d0 <_ZSteqRKSaIcES1_>
   10785:	89 c3                	mov    %eax,%ebx
   10787:	48 8d 45 dd          	lea    -0x23(%rbp),%rax
   1078b:	48 89 c7             	mov    %rax,%rdi
   1078e:	e8 cd 9c ff ff       	call   a460 <_ZNSaIcED1Ev@plt>
   10793:	48 8d 45 de          	lea    -0x22(%rbp),%rax
   10797:	48 89 c7             	mov    %rax,%rdi
   1079a:	e8 c1 9c ff ff       	call   a460 <_ZNSaIcED1Ev@plt>
   1079f:	84 db                	test   %bl,%bl
   107a1:	74 04                	je     107a7 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x9d>
   107a3:	c6 45 df 01          	movb   $0x1,-0x21(%rbp)
   107a7:	80 7d df 00          	cmpb   $0x0,-0x21(%rbp)
   107ab:	0f 84 87 00 00 00    	je     10838 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x12e>
   107b1:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   107b5:	48 89 c7             	mov    %rax,%rdi
   107b8:	e8 d3 99 ff ff       	call   a190 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE4sizeEv@plt>
   107bd:	48 89 c3             	mov    %rax,%rbx
   107c0:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   107c4:	48 89 c7             	mov    %rax,%rdi
   107c7:	e8 c4 99 ff ff       	call   a190 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE4sizeEv@plt>
   107cc:	48 01 d8             	add    %rbx,%rax
   107cf:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   107d3:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   107d7:	48 89 c7             	mov    %rax,%rdi
   107da:	e8 d1 95 ff ff       	call   9db0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE8capacityEv@plt>
   107df:	48 39 45 e0          	cmp    %rax,-0x20(%rbp)
   107e3:	76 19                	jbe    107fe <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0xf4>
   107e5:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   107e9:	48 89 c7             	mov    %rax,%rdi
   107ec:	e8 bf 95 ff ff       	call   9db0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE8capacityEv@plt>
   107f1:	48 39 45 e0          	cmp    %rax,-0x20(%rbp)
   107f5:	77 07                	ja     107fe <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0xf4>
   107f7:	b8 01 00 00 00       	mov    $0x1,%eax
   107fc:	eb 05                	jmp    10803 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0xf9>
   107fe:	b8 00 00 00 00       	mov    $0x0,%eax
   10803:	84 c0                	test   %al,%al
   10805:	74 31                	je     10838 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x12e>
   10807:	48 8b 55 c0          	mov    -0x40(%rbp),%rdx
   1080b:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   1080f:	be 00 00 00 00       	mov    $0x0,%esi
   10814:	48 89 c7             	mov    %rax,%rdi
   10817:	e8 14 a0 ff ff       	call   a830 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmRKS4_@plt>
   1081c:	48 89 c7             	mov    %rax,%rdi
   1081f:	e8 c3 00 00 00       	call   108e7 <_ZSt4moveIRNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEONSt16remove_referenceIT_E4typeEOS8_>
   10824:	48 89 c2             	mov    %rax,%rdx
   10827:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   1082b:	48 89 d6             	mov    %rdx,%rsi
   1082e:	48 89 c7             	mov    %rax,%rdi
   10831:	e8 9a 98 ff ff       	call   a0d0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@plt>
   10836:	eb 2d                	jmp    10865 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x15b>
   10838:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
   1083c:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   10840:	48 89 d6             	mov    %rdx,%rsi
   10843:	48 89 c7             	mov    %rax,%rdi
   10846:	e8 75 a2 ff ff       	call   aac0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendERKS4_@plt>
   1084b:	48 89 c7             	mov    %rax,%rdi
   1084e:	e8 94 00 00 00       	call   108e7 <_ZSt4moveIRNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEONSt16remove_referenceIT_E4typeEOS8_>
   10853:	48 89 c2             	mov    %rax,%rdx
   10856:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   1085a:	48 89 d6             	mov    %rdx,%rsi
   1085d:	48 89 c7             	mov    %rax,%rdi
   10860:	e8 6b 98 ff ff       	call   a0d0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@plt>
   10865:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10869:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   10870:	00 00 
   10872:	74 05                	je     10879 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_+0x16f>
   10874:	e8 57 9c ff ff       	call   a4d0 <__stack_chk_fail@plt>
   10879:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   1087d:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   10881:	c9                   	leave  
   10882:	c3                   	ret    

0000000000010883 <_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEPKS5_OS8_>:
   10883:	f3 0f 1e fa          	endbr64 
   10887:	55                   	push   %rbp
   10888:	48 89 e5             	mov    %rsp,%rbp
   1088b:	48 83 ec 20          	sub    $0x20,%rsp
   1088f:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10893:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   10897:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
   1089b:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   1089f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   108a3:	be 00 00 00 00       	mov    $0x0,%esi
   108a8:	48 89 c7             	mov    %rax,%rdi
   108ab:	e8 30 96 ff ff       	call   9ee0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmPKc@plt>
   108b0:	48 89 c7             	mov    %rax,%rdi
   108b3:	e8 2f 00 00 00       	call   108e7 <_ZSt4moveIRNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEONSt16remove_referenceIT_E4typeEOS8_>
   108b8:	48 89 c2             	mov    %rax,%rdx
   108bb:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   108bf:	48 89 d6             	mov    %rdx,%rsi
   108c2:	48 89 c7             	mov    %rax,%rdi
   108c5:	e8 06 98 ff ff       	call   a0d0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@plt>
   108ca:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   108ce:	c9                   	leave  
   108cf:	c3                   	ret    

00000000000108d0 <_ZSteqRKSaIcES1_>:
   108d0:	f3 0f 1e fa          	endbr64 
   108d4:	55                   	push   %rbp
   108d5:	48 89 e5             	mov    %rsp,%rbp
   108d8:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   108dc:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   108e0:	b8 01 00 00 00       	mov    $0x1,%eax
   108e5:	5d                   	pop    %rbp
   108e6:	c3                   	ret    

00000000000108e7 <_ZSt4moveIRNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEONSt16remove_referenceIT_E4typeEOS8_>:
   108e7:	f3 0f 1e fa          	endbr64 
   108eb:	55                   	push   %rbp
   108ec:	48 89 e5             	mov    %rsp,%rbp
   108ef:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   108f3:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   108f7:	5d                   	pop    %rbp
   108f8:	c3                   	ret    
   108f9:	90                   	nop

00000000000108fa <_ZN3MPI9IntercommD1Ev>:
   108fa:	f3 0f 1e fa          	endbr64 
   108fe:	55                   	push   %rbp
   108ff:	48 89 e5             	mov    %rsp,%rbp
   10902:	48 83 ec 10          	sub    $0x10,%rsp
   10906:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   1090a:	48 8d 15 8f ae 00 00 	lea    0xae8f(%rip),%rdx        # 1b7a0 <_ZTVN3MPI9IntercommE+0x10>
   10911:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10915:	48 89 10             	mov    %rdx,(%rax)
   10918:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1091c:	48 89 c7             	mov    %rax,%rdi
   1091f:	e8 c6 b3 ff ff       	call   bcea <_ZN3MPI4CommD1Ev>
   10924:	90                   	nop
   10925:	c9                   	leave  
   10926:	c3                   	ret    
   10927:	90                   	nop

0000000000010928 <_ZN3MPI9IntercommD0Ev>:
   10928:	f3 0f 1e fa          	endbr64 
   1092c:	55                   	push   %rbp
   1092d:	48 89 e5             	mov    %rsp,%rbp
   10930:	48 83 ec 10          	sub    $0x10,%rsp
   10934:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10938:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1093c:	48 89 c7             	mov    %rax,%rdi
   1093f:	e8 b6 ff ff ff       	call   108fa <_ZN3MPI9IntercommD1Ev>
   10944:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10948:	be 10 00 00 00       	mov    $0x10,%esi
   1094d:	48 89 c7             	mov    %rax,%rdi
   10950:	e8 9b 9a ff ff       	call   a3f0 <_ZdlPvm@plt>
   10955:	c9                   	leave  
   10956:	c3                   	ret    
   10957:	90                   	nop

0000000000010958 <_ZN3MPI9GraphcommD1Ev>:
   10958:	f3 0f 1e fa          	endbr64 
   1095c:	55                   	push   %rbp
   1095d:	48 89 e5             	mov    %rsp,%rbp
   10960:	48 83 ec 10          	sub    $0x10,%rsp
   10964:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10968:	48 8d 15 29 b0 00 00 	lea    0xb029(%rip),%rdx        # 1b998 <_ZTVN3MPI9GraphcommE+0x10>
   1096f:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10973:	48 89 10             	mov    %rdx,(%rax)
   10976:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1097a:	48 89 c7             	mov    %rax,%rdi
   1097d:	e8 f4 b3 ff ff       	call   bd76 <_ZN3MPI9IntracommD1Ev>
   10982:	90                   	nop
   10983:	c9                   	leave  
   10984:	c3                   	ret    
   10985:	90                   	nop

0000000000010986 <_ZN3MPI9GraphcommD0Ev>:
   10986:	f3 0f 1e fa          	endbr64 
   1098a:	55                   	push   %rbp
   1098b:	48 89 e5             	mov    %rsp,%rbp
   1098e:	48 83 ec 10          	sub    $0x10,%rsp
   10992:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10996:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1099a:	48 89 c7             	mov    %rax,%rdi
   1099d:	e8 b6 ff ff ff       	call   10958 <_ZN3MPI9GraphcommD1Ev>
   109a2:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   109a6:	be 10 00 00 00       	mov    $0x10,%esi
   109ab:	48 89 c7             	mov    %rax,%rdi
   109ae:	e8 3d 9a ff ff       	call   a3f0 <_ZdlPvm@plt>
   109b3:	c9                   	leave  
   109b4:	c3                   	ret    
   109b5:	90                   	nop

00000000000109b6 <_ZN3MPI8CartcommD1Ev>:
   109b6:	f3 0f 1e fa          	endbr64 
   109ba:	55                   	push   %rbp
   109bb:	48 89 e5             	mov    %rsp,%rbp
   109be:	48 83 ec 10          	sub    $0x10,%rsp
   109c2:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   109c6:	48 8d 15 2b b2 00 00 	lea    0xb22b(%rip),%rdx        # 1bbf8 <_ZTVN3MPI8CartcommE+0x10>
   109cd:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   109d1:	48 89 10             	mov    %rdx,(%rax)
   109d4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   109d8:	48 89 c7             	mov    %rax,%rdi
   109db:	e8 96 b3 ff ff       	call   bd76 <_ZN3MPI9IntracommD1Ev>
   109e0:	90                   	nop
   109e1:	c9                   	leave  
   109e2:	c3                   	ret    
   109e3:	90                   	nop

00000000000109e4 <_ZN3MPI8CartcommD0Ev>:
   109e4:	f3 0f 1e fa          	endbr64 
   109e8:	55                   	push   %rbp
   109e9:	48 89 e5             	mov    %rsp,%rbp
   109ec:	48 83 ec 10          	sub    $0x10,%rsp
   109f0:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   109f4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   109f8:	48 89 c7             	mov    %rax,%rdi
   109fb:	e8 b6 ff ff ff       	call   109b6 <_ZN3MPI8CartcommD1Ev>
   10a00:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10a04:	be 10 00 00 00       	mov    $0x10,%esi
   10a09:	48 89 c7             	mov    %rax,%rdi
   10a0c:	e8 df 99 ff ff       	call   a3f0 <_ZdlPvm@plt>
   10a11:	c9                   	leave  
   10a12:	c3                   	ret    

0000000000010a13 <_Z14updateVelocityP10particle_td>:
   10a13:	f3 0f 1e fa          	endbr64 
   10a17:	55                   	push   %rbp
   10a18:	48 89 e5             	mov    %rsp,%rbp
   10a1b:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   10a1f:	f2 0f 11 45 d0       	movsd  %xmm0,-0x30(%rbp)
   10a24:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10a28:	0f b6 40 38          	movzbl 0x38(%rax),%eax
   10a2c:	83 f0 01             	xor    $0x1,%eax
   10a2f:	84 c0                	test   %al,%al
   10a31:	75 64                	jne    10a97 <_Z14updateVelocityP10particle_td+0x84>
   10a33:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10a37:	48 8b 40 18          	mov    0x18(%rax),%rax
   10a3b:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   10a3f:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10a43:	48 8b 40 10          	mov    0x10(%rax),%rax
   10a47:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   10a4b:	c7 45 ec 00 00 00 00 	movl   $0x0,-0x14(%rbp)
   10a52:	eb 35                	jmp    10a89 <_Z14updateVelocityP10particle_td+0x76>
   10a54:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   10a58:	48 8d 50 08          	lea    0x8(%rax),%rdx
   10a5c:	48 89 55 f0          	mov    %rdx,-0x10(%rbp)
   10a60:	f2 0f 10 00          	movsd  (%rax),%xmm0
   10a64:	66 0f 28 c8          	movapd %xmm0,%xmm1
   10a68:	f2 0f 59 4d d0       	mulsd  -0x30(%rbp),%xmm1
   10a6d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10a71:	48 8d 50 08          	lea    0x8(%rax),%rdx
   10a75:	48 89 55 f8          	mov    %rdx,-0x8(%rbp)
   10a79:	f2 0f 10 00          	movsd  (%rax),%xmm0
   10a7d:	f2 0f 58 c1          	addsd  %xmm1,%xmm0
   10a81:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   10a85:	83 45 ec 01          	addl   $0x1,-0x14(%rbp)
   10a89:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10a8d:	8b 40 30             	mov    0x30(%rax),%eax
   10a90:	39 45 ec             	cmp    %eax,-0x14(%rbp)
   10a93:	7c bf                	jl     10a54 <_Z14updateVelocityP10particle_td+0x41>
   10a95:	eb 01                	jmp    10a98 <_Z14updateVelocityP10particle_td+0x85>
   10a97:	90                   	nop
   10a98:	5d                   	pop    %rbp
   10a99:	c3                   	ret    

0000000000010a9a <_Z14updatePositionP10particle_td>:
   10a9a:	f3 0f 1e fa          	endbr64 
   10a9e:	55                   	push   %rbp
   10a9f:	48 89 e5             	mov    %rsp,%rbp
   10aa2:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   10aa6:	f2 0f 11 45 d0       	movsd  %xmm0,-0x30(%rbp)
   10aab:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10aaf:	0f b6 40 38          	movzbl 0x38(%rax),%eax
   10ab3:	83 f0 01             	xor    $0x1,%eax
   10ab6:	84 c0                	test   %al,%al
   10ab8:	75 64                	jne    10b1e <_Z14updatePositionP10particle_td+0x84>
   10aba:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10abe:	48 8b 40 08          	mov    0x8(%rax),%rax
   10ac2:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   10ac6:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10aca:	48 8b 40 10          	mov    0x10(%rax),%rax
   10ace:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   10ad2:	c7 45 ec 00 00 00 00 	movl   $0x0,-0x14(%rbp)
   10ad9:	eb 35                	jmp    10b10 <_Z14updatePositionP10particle_td+0x76>
   10adb:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10adf:	48 8d 50 08          	lea    0x8(%rax),%rdx
   10ae3:	48 89 55 f8          	mov    %rdx,-0x8(%rbp)
   10ae7:	f2 0f 10 00          	movsd  (%rax),%xmm0
   10aeb:	66 0f 28 c8          	movapd %xmm0,%xmm1
   10aef:	f2 0f 59 4d d0       	mulsd  -0x30(%rbp),%xmm1
   10af4:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   10af8:	48 8d 50 08          	lea    0x8(%rax),%rdx
   10afc:	48 89 55 f0          	mov    %rdx,-0x10(%rbp)
   10b00:	f2 0f 10 00          	movsd  (%rax),%xmm0
   10b04:	f2 0f 58 c1          	addsd  %xmm1,%xmm0
   10b08:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   10b0c:	83 45 ec 01          	addl   $0x1,-0x14(%rbp)
   10b10:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   10b14:	8b 40 30             	mov    0x30(%rax),%eax
   10b17:	39 45 ec             	cmp    %eax,-0x14(%rbp)
   10b1a:	7c bf                	jl     10adb <_Z14updatePositionP10particle_td+0x41>
   10b1c:	eb 01                	jmp    10b1f <_Z14updatePositionP10particle_td+0x85>
   10b1e:	90                   	nop
   10b1f:	5d                   	pop    %rbp
   10b20:	c3                   	ret    

0000000000010b21 <_Z18updateAccelerationP10particle_tS0_>:
   10b21:	f3 0f 1e fa          	endbr64 
   10b25:	55                   	push   %rbp
   10b26:	48 89 e5             	mov    %rsp,%rbp
   10b29:	48 83 ec 20          	sub    $0x20,%rsp
   10b2d:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   10b31:	48 89 75 e0          	mov    %rsi,-0x20(%rbp)
   10b35:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10b39:	0f b6 40 38          	movzbl 0x38(%rax),%eax
   10b3d:	83 f0 01             	xor    $0x1,%eax
   10b40:	84 c0                	test   %al,%al
   10b42:	74 3a                	je     10b7e <_Z18updateAccelerationP10particle_tS0_+0x5d>
   10b44:	c7 45 f8 00 00 00 00 	movl   $0x0,-0x8(%rbp)
   10b4b:	eb 20                	jmp    10b6d <_Z18updateAccelerationP10particle_tS0_+0x4c>
   10b4d:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10b51:	48 8b 50 18          	mov    0x18(%rax),%rdx
   10b55:	8b 45 f8             	mov    -0x8(%rbp),%eax
   10b58:	48 98                	cltq   
   10b5a:	48 c1 e0 03          	shl    $0x3,%rax
   10b5e:	48 01 d0             	add    %rdx,%rax
   10b61:	66 0f ef c0          	pxor   %xmm0,%xmm0
   10b65:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   10b69:	83 45 f8 01          	addl   $0x1,-0x8(%rbp)
   10b6d:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10b71:	8b 40 30             	mov    0x30(%rax),%eax
   10b74:	39 45 f8             	cmp    %eax,-0x8(%rbp)
   10b77:	7c d4                	jl     10b4d <_Z18updateAccelerationP10particle_tS0_+0x2c>
   10b79:	e9 a8 00 00 00       	jmp    10c26 <_Z18updateAccelerationP10particle_tS0_+0x105>
   10b7e:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   10b82:	0f b6 40 38          	movzbl 0x38(%rax),%eax
   10b86:	83 f0 01             	xor    $0x1,%eax
   10b89:	84 c0                	test   %al,%al
   10b8b:	74 37                	je     10bc4 <_Z18updateAccelerationP10particle_tS0_+0xa3>
   10b8d:	c7 45 fc 00 00 00 00 	movl   $0x0,-0x4(%rbp)
   10b94:	eb 20                	jmp    10bb6 <_Z18updateAccelerationP10particle_tS0_+0x95>
   10b96:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   10b9a:	48 8b 50 18          	mov    0x18(%rax),%rdx
   10b9e:	8b 45 fc             	mov    -0x4(%rbp),%eax
   10ba1:	48 98                	cltq   
   10ba3:	48 c1 e0 03          	shl    $0x3,%rax
   10ba7:	48 01 d0             	add    %rdx,%rax
   10baa:	66 0f ef c0          	pxor   %xmm0,%xmm0
   10bae:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   10bb2:	83 45 fc 01          	addl   $0x1,-0x4(%rbp)
   10bb6:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   10bba:	8b 40 30             	mov    0x30(%rax),%eax
   10bbd:	39 45 fc             	cmp    %eax,-0x4(%rbp)
   10bc0:	7c d4                	jl     10b96 <_Z18updateAccelerationP10particle_tS0_+0x75>
   10bc2:	eb 62                	jmp    10c26 <_Z18updateAccelerationP10particle_tS0_+0x105>
   10bc4:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10bc8:	48 8b 50 28          	mov    0x28(%rax),%rdx
   10bcc:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   10bd0:	48 8b 40 28          	mov    0x28(%rax),%rax
   10bd4:	48 39 c2             	cmp    %rax,%rdx
   10bd7:	74 28                	je     10c01 <_Z18updateAccelerationP10particle_tS0_+0xe0>
   10bd9:	48 8d 05 a8 5b 00 00 	lea    0x5ba8(%rip),%rax        # 16788 <_ZStL19piecewise_construct+0x8>
   10be0:	48 89 c1             	mov    %rax,%rcx
   10be3:	ba 2b 00 00 00       	mov    $0x2b,%edx
   10be8:	48 8d 05 d1 5b 00 00 	lea    0x5bd1(%rip),%rax        # 167c0 <_ZStL19piecewise_construct+0x40>
   10bef:	48 89 c6             	mov    %rax,%rsi
   10bf2:	48 8d 05 1f 5c 00 00 	lea    0x5c1f(%rip),%rax        # 16818 <_ZStL19piecewise_construct+0x98>
   10bf9:	48 89 c7             	mov    %rax,%rdi
   10bfc:	e8 2f 95 ff ff       	call   a130 <__assert_fail@plt>
   10c01:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10c05:	48 8b 40 28          	mov    0x28(%rax),%rax
   10c09:	48 85 c0             	test   %rax,%rax
   10c0c:	74 18                	je     10c26 <_Z18updateAccelerationP10particle_tS0_+0x105>
   10c0e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10c12:	48 8b 48 28          	mov    0x28(%rax),%rcx
   10c16:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
   10c1a:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10c1e:	48 89 d6             	mov    %rdx,%rsi
   10c21:	48 89 c7             	mov    %rax,%rdi
   10c24:	ff d1                	call   *%rcx
   10c26:	c9                   	leave  
   10c27:	c3                   	ret    

0000000000010c28 <_Z28registerUpdateAccelerationFnP10particle_tPFvS0_S0_E>:
   10c28:	f3 0f 1e fa          	endbr64 
   10c2c:	55                   	push   %rbp
   10c2d:	48 89 e5             	mov    %rsp,%rbp
   10c30:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10c34:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   10c38:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10c3c:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   10c40:	48 89 50 28          	mov    %rdx,0x28(%rax)
   10c44:	90                   	nop
   10c45:	5d                   	pop    %rbp
   10c46:	c3                   	ret    

0000000000010c47 <_Z29registerUpdateAccelerationFnsP10particle_tiPFvS0_S0_E>:
   10c47:	f3 0f 1e fa          	endbr64 
   10c4b:	55                   	push   %rbp
   10c4c:	48 89 e5             	mov    %rsp,%rbp
   10c4f:	48 83 ec 28          	sub    $0x28,%rsp
   10c53:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   10c57:	89 75 e4             	mov    %esi,-0x1c(%rbp)
   10c5a:	48 89 55 d8          	mov    %rdx,-0x28(%rbp)
   10c5e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   10c62:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   10c66:	c7 45 f4 00 00 00 00 	movl   $0x0,-0xc(%rbp)
   10c6d:	eb 1f                	jmp    10c8e <_Z29registerUpdateAccelerationFnsP10particle_tiPFvS0_S0_E+0x47>
   10c6f:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10c73:	48 8d 50 40          	lea    0x40(%rax),%rdx
   10c77:	48 89 55 f8          	mov    %rdx,-0x8(%rbp)
   10c7b:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
   10c7f:	48 89 d6             	mov    %rdx,%rsi
   10c82:	48 89 c7             	mov    %rax,%rdi
   10c85:	e8 9e ff ff ff       	call   10c28 <_Z28registerUpdateAccelerationFnP10particle_tPFvS0_S0_E>
   10c8a:	83 45 f4 01          	addl   $0x1,-0xc(%rbp)
   10c8e:	8b 45 f4             	mov    -0xc(%rbp),%eax
   10c91:	3b 45 e4             	cmp    -0x1c(%rbp),%eax
   10c94:	7c d9                	jl     10c6f <_Z29registerUpdateAccelerationFnsP10particle_tiPFvS0_S0_E+0x28>
   10c96:	90                   	nop
   10c97:	90                   	nop
   10c98:	c9                   	leave  
   10c99:	c3                   	ret    

0000000000010c9a <_Z29registerUpdateAccelerationFnsP17chunk_particles_tPFvP10particle_tS2_E>:
   10c9a:	f3 0f 1e fa          	endbr64 
   10c9e:	55                   	push   %rbp
   10c9f:	48 89 e5             	mov    %rsp,%rbp
   10ca2:	48 83 ec 10          	sub    $0x10,%rsp
   10ca6:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   10caa:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   10cae:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10cb2:	8b 48 08             	mov    0x8(%rax),%ecx
   10cb5:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   10cb9:	48 8b 00             	mov    (%rax),%rax
   10cbc:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   10cc0:	89 ce                	mov    %ecx,%esi
   10cc2:	48 89 c7             	mov    %rax,%rdi
   10cc5:	e8 7d ff ff ff       	call   10c47 <_Z29registerUpdateAccelerationFnsP10particle_tiPFvS0_S0_E>
   10cca:	90                   	nop
   10ccb:	c9                   	leave  
   10ccc:	c3                   	ret    

0000000000010ccd <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi>:
   10ccd:	f3 0f 1e fa          	endbr64 
   10cd1:	55                   	push   %rbp
   10cd2:	48 89 e5             	mov    %rsp,%rbp
   10cd5:	53                   	push   %rbx
   10cd6:	48 81 ec 48 04 00 00 	sub    $0x448,%rsp
   10cdd:	48 89 bd c8 fb ff ff 	mov    %rdi,-0x438(%rbp)
   10ce4:	48 89 b5 c0 fb ff ff 	mov    %rsi,-0x440(%rbp)
   10ceb:	48 89 95 b8 fb ff ff 	mov    %rdx,-0x448(%rbp)
   10cf2:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   10cf9:	00 00 
   10cfb:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   10cff:	31 c0                	xor    %eax,%eax
   10d01:	48 8b 8d c8 fb ff ff 	mov    -0x438(%rbp),%rcx
   10d08:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   10d0f:	ba 08 00 00 00       	mov    $0x8,%edx
   10d14:	48 89 ce             	mov    %rcx,%rsi
   10d17:	48 89 c7             	mov    %rax,%rdi
   10d1a:	e8 51 9b ff ff       	call   a870 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEEC1ERKNSt7__cxx1112basic_stringIcS1_SaIcEEESt13_Ios_Openmode@plt>
   10d1f:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   10d26:	48 89 c7             	mov    %rax,%rdi
   10d29:	e8 22 9b ff ff       	call   a850 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEE7is_openEv@plt>
   10d2e:	83 f0 01             	xor    $0x1,%eax
   10d31:	84 c0                	test   %al,%al
   10d33:	74 63                	je     10d98 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0xcb>
   10d35:	48 8b 95 c8 fb ff ff 	mov    -0x438(%rbp),%rdx
   10d3c:	48 8d 85 50 fc ff ff 	lea    -0x3b0(%rbp),%rax
   10d43:	48 89 d6             	mov    %rdx,%rsi
   10d46:	48 89 c7             	mov    %rax,%rdi
   10d49:	e8 72 92 ff ff       	call   9fc0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
   10d4e:	48 8b 05 cb d4 00 00 	mov    0xd4cb(%rip),%rax        # 1e220 <stderr@GLIBC_2.2.5>
   10d55:	48 8d 95 50 fc ff ff 	lea    -0x3b0(%rbp),%rdx
   10d5c:	49 89 d0             	mov    %rdx,%r8
   10d5f:	b9 45 00 00 00       	mov    $0x45,%ecx
   10d64:	48 8d 15 55 5a 00 00 	lea    0x5a55(%rip),%rdx        # 167c0 <_ZStL19piecewise_construct+0x40>
   10d6b:	48 8d 35 de 5a 00 00 	lea    0x5ade(%rip),%rsi        # 16850 <_ZStL19piecewise_construct+0xd0>
   10d72:	48 89 c7             	mov    %rax,%rdi
   10d75:	b8 00 00 00 00       	mov    $0x0,%eax
   10d7a:	e8 d1 99 ff ff       	call   a750 <fprintf@plt>
   10d7f:	48 8d 85 50 fc ff ff 	lea    -0x3b0(%rbp),%rax
   10d86:	48 89 c7             	mov    %rax,%rdi
   10d89:	e8 72 93 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   10d8e:	bf 01 00 00 00       	mov    $0x1,%edi
   10d93:	e8 98 97 ff ff       	call   a530 <exit@plt>
   10d98:	c7 85 d0 fb ff ff ff 	movl   $0xffffffff,-0x430(%rbp)
   10d9f:	ff ff ff 
   10da2:	48 8d 85 30 fc ff ff 	lea    -0x3d0(%rbp),%rax
   10da9:	48 89 c7             	mov    %rax,%rdi
   10dac:	e8 df 98 ff ff       	call   a690 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1Ev@plt>
   10db1:	48 8d 95 30 fc ff ff 	lea    -0x3d0(%rbp),%rdx
   10db8:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   10dbf:	48 89 d6             	mov    %rdx,%rsi
   10dc2:	48 89 c7             	mov    %rax,%rdi
   10dc5:	e8 96 90 ff ff       	call   9e60 <_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE@plt>
   10dca:	48 8d 85 30 fc ff ff 	lea    -0x3d0(%rbp),%rax
   10dd1:	48 89 c7             	mov    %rax,%rdi
   10dd4:	e8 c7 91 ff ff       	call   9fa0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE5c_strEv@plt>
   10dd9:	48 89 c7             	mov    %rax,%rdi
   10ddc:	48 8d 95 d0 fb ff ff 	lea    -0x430(%rbp),%rdx
   10de3:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   10dea:	48 89 d1             	mov    %rdx,%rcx
   10ded:	48 89 c2             	mov    %rax,%rdx
   10df0:	48 8d 05 7c 5a 00 00 	lea    0x5a7c(%rip),%rax        # 16873 <_ZStL19piecewise_construct+0xf3>
   10df7:	48 89 c6             	mov    %rax,%rsi
   10dfa:	b8 00 00 00 00       	mov    $0x0,%eax
   10dff:	e8 ac 95 ff ff       	call   a3b0 <__isoc99_sscanf@plt>
   10e04:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   10e0b:	8b 00                	mov    (%rax),%eax
   10e0d:	48 98                	cltq   
   10e0f:	48 c1 e0 06          	shl    $0x6,%rax
   10e13:	48 89 c7             	mov    %rax,%rdi
   10e16:	e8 05 98 ff ff       	call   a620 <malloc@plt>
   10e1b:	48 89 85 00 fc ff ff 	mov    %rax,-0x400(%rbp)
   10e22:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   10e29:	8b 10                	mov    (%rax),%edx
   10e2b:	8b 85 d0 fb ff ff    	mov    -0x430(%rbp),%eax
   10e31:	0f af c2             	imul   %edx,%eax
   10e34:	48 98                	cltq   
   10e36:	48 c1 e0 03          	shl    $0x3,%rax
   10e3a:	48 89 c7             	mov    %rax,%rdi
   10e3d:	e8 de 97 ff ff       	call   a620 <malloc@plt>
   10e42:	48 89 85 08 fc ff ff 	mov    %rax,-0x3f8(%rbp)
   10e49:	48 83 bd 08 fc ff ff 	cmpq   $0x0,-0x3f8(%rbp)
   10e50:	00 
   10e51:	75 28                	jne    10e7b <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x1ae>
   10e53:	48 8d 05 26 5a 00 00 	lea    0x5a26(%rip),%rax        # 16880 <_ZStL19piecewise_construct+0x100>
   10e5a:	48 89 c1             	mov    %rax,%rcx
   10e5d:	ba 4e 00 00 00       	mov    $0x4e,%edx
   10e62:	48 8d 05 57 59 00 00 	lea    0x5957(%rip),%rax        # 167c0 <_ZStL19piecewise_construct+0x40>
   10e69:	48 89 c6             	mov    %rax,%rsi
   10e6c:	48 8d 05 42 5a 00 00 	lea    0x5a42(%rip),%rax        # 168b5 <_ZStL19piecewise_construct+0x135>
   10e73:	48 89 c7             	mov    %rax,%rdi
   10e76:	e8 b5 92 ff ff       	call   a130 <__assert_fail@plt>
   10e7b:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   10e82:	8b 10                	mov    (%rax),%edx
   10e84:	8b 85 d0 fb ff ff    	mov    -0x430(%rbp),%eax
   10e8a:	0f af c2             	imul   %edx,%eax
   10e8d:	48 98                	cltq   
   10e8f:	48 c1 e0 03          	shl    $0x3,%rax
   10e93:	48 89 c7             	mov    %rax,%rdi
   10e96:	e8 85 97 ff ff       	call   a620 <malloc@plt>
   10e9b:	48 89 85 10 fc ff ff 	mov    %rax,-0x3f0(%rbp)
   10ea2:	48 83 bd 10 fc ff ff 	cmpq   $0x0,-0x3f0(%rbp)
   10ea9:	00 
   10eaa:	75 28                	jne    10ed4 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x207>
   10eac:	48 8d 05 cd 59 00 00 	lea    0x59cd(%rip),%rax        # 16880 <_ZStL19piecewise_construct+0x100>
   10eb3:	48 89 c1             	mov    %rax,%rcx
   10eb6:	ba 50 00 00 00       	mov    $0x50,%edx
   10ebb:	48 8d 05 fe 58 00 00 	lea    0x58fe(%rip),%rax        # 167c0 <_ZStL19piecewise_construct+0x40>
   10ec2:	48 89 c6             	mov    %rax,%rsi
   10ec5:	48 8d 05 f8 59 00 00 	lea    0x59f8(%rip),%rax        # 168c4 <_ZStL19piecewise_construct+0x144>
   10ecc:	48 89 c7             	mov    %rax,%rdi
   10ecf:	e8 5c 92 ff ff       	call   a130 <__assert_fail@plt>
   10ed4:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   10edb:	8b 10                	mov    (%rax),%edx
   10edd:	8b 85 d0 fb ff ff    	mov    -0x430(%rbp),%eax
   10ee3:	0f af c2             	imul   %edx,%eax
   10ee6:	48 98                	cltq   
   10ee8:	48 c1 e0 03          	shl    $0x3,%rax
   10eec:	48 89 c7             	mov    %rax,%rdi
   10eef:	e8 2c 97 ff ff       	call   a620 <malloc@plt>
   10ef4:	48 89 85 18 fc ff ff 	mov    %rax,-0x3e8(%rbp)
   10efb:	48 83 bd 18 fc ff ff 	cmpq   $0x0,-0x3e8(%rbp)
   10f02:	00 
   10f03:	75 28                	jne    10f2d <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x260>
   10f05:	48 8d 05 74 59 00 00 	lea    0x5974(%rip),%rax        # 16880 <_ZStL19piecewise_construct+0x100>
   10f0c:	48 89 c1             	mov    %rax,%rcx
   10f0f:	ba 52 00 00 00       	mov    $0x52,%edx
   10f14:	48 8d 05 a5 58 00 00 	lea    0x58a5(%rip),%rax        # 167c0 <_ZStL19piecewise_construct+0x40>
   10f1b:	48 89 c6             	mov    %rax,%rsi
   10f1e:	48 8d 05 ae 59 00 00 	lea    0x59ae(%rip),%rax        # 168d3 <_ZStL19piecewise_construct+0x153>
   10f25:	48 89 c7             	mov    %rax,%rdi
   10f28:	e8 03 92 ff ff       	call   a130 <__assert_fail@plt>
   10f2d:	48 c7 85 f0 fb ff ff 	movq   $0x0,-0x410(%rbp)
   10f34:	00 00 00 00 
   10f38:	48 8b 85 00 fc ff ff 	mov    -0x400(%rbp),%rax
   10f3f:	48 89 85 f8 fb ff ff 	mov    %rax,-0x408(%rbp)
   10f46:	c7 85 d4 fb ff ff 00 	movl   $0x0,-0x42c(%rbp)
   10f4d:	00 00 00 
   10f50:	e9 da 00 00 00       	jmp    1102f <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x362>
   10f55:	8b 85 d0 fb ff ff    	mov    -0x430(%rbp),%eax
   10f5b:	0f af 85 d4 fb ff ff 	imul   -0x42c(%rbp),%eax
   10f62:	48 98                	cltq   
   10f64:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   10f6b:	00 
   10f6c:	48 8b 85 08 fc ff ff 	mov    -0x3f8(%rbp),%rax
   10f73:	48 01 c2             	add    %rax,%rdx
   10f76:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   10f7d:	48 89 50 08          	mov    %rdx,0x8(%rax)
   10f81:	8b 85 d0 fb ff ff    	mov    -0x430(%rbp),%eax
   10f87:	0f af 85 d4 fb ff ff 	imul   -0x42c(%rbp),%eax
   10f8e:	48 98                	cltq   
   10f90:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   10f97:	00 
   10f98:	48 8b 85 10 fc ff ff 	mov    -0x3f0(%rbp),%rax
   10f9f:	48 01 c2             	add    %rax,%rdx
   10fa2:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   10fa9:	48 89 50 10          	mov    %rdx,0x10(%rax)
   10fad:	8b 85 d0 fb ff ff    	mov    -0x430(%rbp),%eax
   10fb3:	0f af 85 d4 fb ff ff 	imul   -0x42c(%rbp),%eax
   10fba:	48 98                	cltq   
   10fbc:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   10fc3:	00 
   10fc4:	48 8b 85 18 fc ff ff 	mov    -0x3e8(%rbp),%rax
   10fcb:	48 01 c2             	add    %rax,%rdx
   10fce:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   10fd5:	48 89 50 18          	mov    %rdx,0x18(%rax)
   10fd9:	8b 95 d0 fb ff ff    	mov    -0x430(%rbp),%edx
   10fdf:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   10fe6:	89 50 30             	mov    %edx,0x30(%rax)
   10fe9:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   10ff0:	48 c7 40 28 00 00 00 	movq   $0x0,0x28(%rax)
   10ff7:	00 
   10ff8:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   10fff:	48 c7 40 20 00 00 00 	movq   $0x0,0x20(%rax)
   11006:	00 
   11007:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   1100e:	c7 40 34 ff ff ff ff 	movl   $0xffffffff,0x34(%rax)
   11015:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   1101c:	c6 40 38 01          	movb   $0x1,0x38(%rax)
   11020:	48 83 85 f8 fb ff ff 	addq   $0x40,-0x408(%rbp)
   11027:	40 
   11028:	83 85 d4 fb ff ff 01 	addl   $0x1,-0x42c(%rbp)
   1102f:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   11036:	8b 00                	mov    (%rax),%eax
   11038:	39 85 d4 fb ff ff    	cmp    %eax,-0x42c(%rbp)
   1103e:	0f 8c 11 ff ff ff    	jl     10f55 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x288>
   11044:	48 8b 85 00 fc ff ff 	mov    -0x400(%rbp),%rax
   1104b:	48 89 85 f8 fb ff ff 	mov    %rax,-0x408(%rbp)
   11052:	48 8d 85 50 fc ff ff 	lea    -0x3b0(%rbp),%rax
   11059:	48 89 c7             	mov    %rax,%rdi
   1105c:	e8 4f 91 ff ff       	call   a1b0 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEC1Ev@plt>
   11061:	c7 85 d8 fb ff ff 00 	movl   $0x0,-0x428(%rbp)
   11068:	00 00 00 
   1106b:	c7 85 dc fb ff ff ff 	movl   $0xffffffff,-0x424(%rbp)
   11072:	ff ff ff 
   11075:	e9 60 03 00 00       	jmp    113da <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x70d>
   1107a:	48 8d 85 30 fc ff ff 	lea    -0x3d0(%rbp),%rax
   11081:	48 89 c7             	mov    %rax,%rdi
   11084:	e8 17 8f ff ff       	call   9fa0 <_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE5c_strEv@plt>
   11089:	48 89 c2             	mov    %rax,%rdx
   1108c:	48 8d 05 4f 58 00 00 	lea    0x584f(%rip),%rax        # 168e2 <_ZStL19piecewise_construct+0x162>
   11093:	48 89 c6             	mov    %rax,%rsi
   11096:	48 89 d7             	mov    %rdx,%rdi
   11099:	e8 72 97 ff ff       	call   a810 <strtok@plt>
   1109e:	48 89 85 20 fc ff ff 	mov    %rax,-0x3e0(%rbp)
   110a5:	48 8b 85 20 fc ff ff 	mov    -0x3e0(%rbp),%rax
   110ac:	48 89 c7             	mov    %rax,%rdi
   110af:	e8 bc 94 ff ff       	call   a570 <atof@plt>
   110b4:	66 48 0f 7e c0       	movq   %xmm0,%rax
   110b9:	48 89 85 28 fc ff ff 	mov    %rax,-0x3d8(%rbp)
   110c0:	f2 0f 10 85 28 fc ff 	movsd  -0x3d8(%rbp),%xmm0
   110c7:	ff 
   110c8:	66 0f 2f 05 28 59 00 	comisd 0x5928(%rip),%xmm0        # 169f8 <_ZStL19piecewise_construct+0x278>
   110cf:	00 
   110d0:	73 0f                	jae    110e1 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x414>
   110d2:	f2 0f 10 85 28 fc ff 	movsd  -0x3d8(%rbp),%xmm0
   110d9:	ff 
   110da:	f2 48 0f 2c c0       	cvttsd2si %xmm0,%rax
   110df:	eb 26                	jmp    11107 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x43a>
   110e1:	f2 0f 10 85 28 fc ff 	movsd  -0x3d8(%rbp),%xmm0
   110e8:	ff 
   110e9:	f2 0f 10 0d 07 59 00 	movsd  0x5907(%rip),%xmm1        # 169f8 <_ZStL19piecewise_construct+0x278>
   110f0:	00 
   110f1:	f2 0f 5c c1          	subsd  %xmm1,%xmm0
   110f5:	f2 48 0f 2c c0       	cvttsd2si %xmm0,%rax
   110fa:	48 ba 00 00 00 00 00 	movabs $0x8000000000000000,%rdx
   11101:	00 00 80 
   11104:	48 31 d0             	xor    %rdx,%rax
   11107:	48 8b 95 f8 fb ff ff 	mov    -0x408(%rbp),%rdx
   1110e:	48 89 02             	mov    %rax,(%rdx)
   11111:	c7 85 e0 fb ff ff 00 	movl   $0x0,-0x420(%rbp)
   11118:	00 00 00 
   1111b:	eb 54                	jmp    11171 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x4a4>
   1111d:	48 8d 05 be 57 00 00 	lea    0x57be(%rip),%rax        # 168e2 <_ZStL19piecewise_construct+0x162>
   11124:	48 89 c6             	mov    %rax,%rsi
   11127:	bf 00 00 00 00       	mov    $0x0,%edi
   1112c:	e8 df 96 ff ff       	call   a810 <strtok@plt>
   11131:	48 89 85 20 fc ff ff 	mov    %rax,-0x3e0(%rbp)
   11138:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   1113f:	48 8b 50 08          	mov    0x8(%rax),%rdx
   11143:	8b 85 e0 fb ff ff    	mov    -0x420(%rbp),%eax
   11149:	48 98                	cltq   
   1114b:	48 c1 e0 03          	shl    $0x3,%rax
   1114f:	48 8d 1c 02          	lea    (%rdx,%rax,1),%rbx
   11153:	48 8b 85 20 fc ff ff 	mov    -0x3e0(%rbp),%rax
   1115a:	48 89 c7             	mov    %rax,%rdi
   1115d:	e8 0e 94 ff ff       	call   a570 <atof@plt>
   11162:	66 48 0f 7e c0       	movq   %xmm0,%rax
   11167:	48 89 03             	mov    %rax,(%rbx)
   1116a:	83 85 e0 fb ff ff 01 	addl   $0x1,-0x420(%rbp)
   11171:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   11178:	8b 40 30             	mov    0x30(%rax),%eax
   1117b:	39 85 e0 fb ff ff    	cmp    %eax,-0x420(%rbp)
   11181:	7c 9a                	jl     1111d <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x450>
   11183:	c7 85 e4 fb ff ff 00 	movl   $0x0,-0x41c(%rbp)
   1118a:	00 00 00 
   1118d:	eb 54                	jmp    111e3 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x516>
   1118f:	48 8d 05 4c 57 00 00 	lea    0x574c(%rip),%rax        # 168e2 <_ZStL19piecewise_construct+0x162>
   11196:	48 89 c6             	mov    %rax,%rsi
   11199:	bf 00 00 00 00       	mov    $0x0,%edi
   1119e:	e8 6d 96 ff ff       	call   a810 <strtok@plt>
   111a3:	48 89 85 20 fc ff ff 	mov    %rax,-0x3e0(%rbp)
   111aa:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   111b1:	48 8b 50 10          	mov    0x10(%rax),%rdx
   111b5:	8b 85 e4 fb ff ff    	mov    -0x41c(%rbp),%eax
   111bb:	48 98                	cltq   
   111bd:	48 c1 e0 03          	shl    $0x3,%rax
   111c1:	48 8d 1c 02          	lea    (%rdx,%rax,1),%rbx
   111c5:	48 8b 85 20 fc ff ff 	mov    -0x3e0(%rbp),%rax
   111cc:	48 89 c7             	mov    %rax,%rdi
   111cf:	e8 9c 93 ff ff       	call   a570 <atof@plt>
   111d4:	66 48 0f 7e c0       	movq   %xmm0,%rax
   111d9:	48 89 03             	mov    %rax,(%rbx)
   111dc:	83 85 e4 fb ff ff 01 	addl   $0x1,-0x41c(%rbp)
   111e3:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   111ea:	8b 40 30             	mov    0x30(%rax),%eax
   111ed:	39 85 e4 fb ff ff    	cmp    %eax,-0x41c(%rbp)
   111f3:	7c 9a                	jl     1118f <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x4c2>
   111f5:	c7 85 e8 fb ff ff 00 	movl   $0x0,-0x418(%rbp)
   111fc:	00 00 00 
   111ff:	eb 54                	jmp    11255 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x588>
   11201:	48 8d 05 da 56 00 00 	lea    0x56da(%rip),%rax        # 168e2 <_ZStL19piecewise_construct+0x162>
   11208:	48 89 c6             	mov    %rax,%rsi
   1120b:	bf 00 00 00 00       	mov    $0x0,%edi
   11210:	e8 fb 95 ff ff       	call   a810 <strtok@plt>
   11215:	48 89 85 20 fc ff ff 	mov    %rax,-0x3e0(%rbp)
   1121c:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   11223:	48 8b 50 18          	mov    0x18(%rax),%rdx
   11227:	8b 85 e8 fb ff ff    	mov    -0x418(%rbp),%eax
   1122d:	48 98                	cltq   
   1122f:	48 c1 e0 03          	shl    $0x3,%rax
   11233:	48 8d 1c 02          	lea    (%rdx,%rax,1),%rbx
   11237:	48 8b 85 20 fc ff ff 	mov    -0x3e0(%rbp),%rax
   1123e:	48 89 c7             	mov    %rax,%rdi
   11241:	e8 2a 93 ff ff       	call   a570 <atof@plt>
   11246:	66 48 0f 7e c0       	movq   %xmm0,%rax
   1124b:	48 89 03             	mov    %rax,(%rbx)
   1124e:	83 85 e8 fb ff ff 01 	addl   $0x1,-0x418(%rbp)
   11255:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   1125c:	8b 40 30             	mov    0x30(%rax),%eax
   1125f:	39 85 e8 fb ff ff    	cmp    %eax,-0x418(%rbp)
   11265:	7c 9a                	jl     11201 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x534>
   11267:	48 8d 05 74 56 00 00 	lea    0x5674(%rip),%rax        # 168e2 <_ZStL19piecewise_construct+0x162>
   1126e:	48 89 c6             	mov    %rax,%rsi
   11271:	bf 00 00 00 00       	mov    $0x0,%edi
   11276:	e8 95 95 ff ff       	call   a810 <strtok@plt>
   1127b:	48 89 85 20 fc ff ff 	mov    %rax,-0x3e0(%rbp)
   11282:	48 8b 85 20 fc ff ff 	mov    -0x3e0(%rbp),%rax
   11289:	48 89 c7             	mov    %rax,%rdi
   1128c:	e8 df 92 ff ff       	call   a570 <atof@plt>
   11291:	66 48 0f 7e c0       	movq   %xmm0,%rax
   11296:	48 89 85 28 fc ff ff 	mov    %rax,-0x3d8(%rbp)
   1129d:	f2 0f 10 85 28 fc ff 	movsd  -0x3d8(%rbp),%xmm0
   112a4:	ff 
   112a5:	f2 0f 2c d0          	cvttsd2si %xmm0,%edx
   112a9:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   112b0:	89 50 34             	mov    %edx,0x34(%rax)
   112b3:	83 bd dc fb ff ff ff 	cmpl   $0xffffffff,-0x424(%rbp)
   112ba:	75 37                	jne    112f3 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x626>
   112bc:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   112c3:	8b 40 34             	mov    0x34(%rax),%eax
   112c6:	89 85 dc fb ff ff    	mov    %eax,-0x424(%rbp)
   112cc:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   112d3:	8b 00                	mov    (%rax),%eax
   112d5:	0f af 85 dc fb ff ff 	imul   -0x424(%rbp),%eax
   112dc:	48 98                	cltq   
   112de:	48 c1 e0 03          	shl    $0x3,%rax
   112e2:	48 89 c7             	mov    %rax,%rdi
   112e5:	e8 36 93 ff ff       	call   a620 <malloc@plt>
   112ea:	48 89 85 f0 fb ff ff 	mov    %rax,-0x410(%rbp)
   112f1:	eb 3a                	jmp    1132d <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x660>
   112f3:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   112fa:	8b 40 34             	mov    0x34(%rax),%eax
   112fd:	39 85 dc fb ff ff    	cmp    %eax,-0x424(%rbp)
   11303:	74 28                	je     1132d <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x660>
   11305:	48 8d 05 74 55 00 00 	lea    0x5574(%rip),%rax        # 16880 <_ZStL19piecewise_construct+0x100>
   1130c:	48 89 c1             	mov    %rax,%rcx
   1130f:	ba 7f 00 00 00       	mov    $0x7f,%edx
   11314:	48 8d 05 a5 54 00 00 	lea    0x54a5(%rip),%rax        # 167c0 <_ZStL19piecewise_construct+0x40>
   1131b:	48 89 c6             	mov    %rax,%rsi
   1131e:	48 8d 05 bf 55 00 00 	lea    0x55bf(%rip),%rax        # 168e4 <_ZStL19piecewise_construct+0x164>
   11325:	48 89 c7             	mov    %rax,%rdi
   11328:	e8 03 8e ff ff       	call   a130 <__assert_fail@plt>
   1132d:	8b 85 d8 fb ff ff    	mov    -0x428(%rbp),%eax
   11333:	0f af 85 dc fb ff ff 	imul   -0x424(%rbp),%eax
   1133a:	48 98                	cltq   
   1133c:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   11343:	00 
   11344:	48 8b 85 f0 fb ff ff 	mov    -0x410(%rbp),%rax
   1134b:	48 01 c2             	add    %rax,%rdx
   1134e:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   11355:	48 89 50 20          	mov    %rdx,0x20(%rax)
   11359:	c7 85 ec fb ff ff 00 	movl   $0x0,-0x414(%rbp)
   11360:	00 00 00 
   11363:	eb 54                	jmp    113b9 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x6ec>
   11365:	48 8d 05 76 55 00 00 	lea    0x5576(%rip),%rax        # 168e2 <_ZStL19piecewise_construct+0x162>
   1136c:	48 89 c6             	mov    %rax,%rsi
   1136f:	bf 00 00 00 00       	mov    $0x0,%edi
   11374:	e8 97 94 ff ff       	call   a810 <strtok@plt>
   11379:	48 89 85 20 fc ff ff 	mov    %rax,-0x3e0(%rbp)
   11380:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   11387:	48 8b 50 20          	mov    0x20(%rax),%rdx
   1138b:	8b 85 ec fb ff ff    	mov    -0x414(%rbp),%eax
   11391:	48 98                	cltq   
   11393:	48 c1 e0 03          	shl    $0x3,%rax
   11397:	48 8d 1c 02          	lea    (%rdx,%rax,1),%rbx
   1139b:	48 8b 85 20 fc ff ff 	mov    -0x3e0(%rbp),%rax
   113a2:	48 89 c7             	mov    %rax,%rdi
   113a5:	e8 c6 91 ff ff       	call   a570 <atof@plt>
   113aa:	66 48 0f 7e c0       	movq   %xmm0,%rax
   113af:	48 89 03             	mov    %rax,(%rbx)
   113b2:	83 85 ec fb ff ff 01 	addl   $0x1,-0x414(%rbp)
   113b9:	48 8b 85 f8 fb ff ff 	mov    -0x408(%rbp),%rax
   113c0:	8b 40 34             	mov    0x34(%rax),%eax
   113c3:	39 85 ec fb ff ff    	cmp    %eax,-0x414(%rbp)
   113c9:	7c 9a                	jl     11365 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x698>
   113cb:	83 85 d8 fb ff ff 01 	addl   $0x1,-0x428(%rbp)
   113d2:	48 83 85 f8 fb ff ff 	addq   $0x40,-0x408(%rbp)
   113d9:	40 
   113da:	48 8d 95 30 fc ff ff 	lea    -0x3d0(%rbp),%rdx
   113e1:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   113e8:	48 89 d6             	mov    %rdx,%rsi
   113eb:	48 89 c7             	mov    %rax,%rdi
   113ee:	e8 6d 8a ff ff       	call   9e60 <_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE@plt>
   113f3:	48 8b 10             	mov    (%rax),%rdx
   113f6:	48 83 ea 18          	sub    $0x18,%rdx
   113fa:	48 8b 12             	mov    (%rdx),%rdx
   113fd:	48 01 d0             	add    %rdx,%rax
   11400:	48 89 c7             	mov    %rax,%rdi
   11403:	e8 88 91 ff ff       	call   a590 <_ZNKSt9basic_iosIcSt11char_traitsIcEEcvbEv@plt>
   11408:	84 c0                	test   %al,%al
   1140a:	74 18                	je     11424 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x757>
   1140c:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   11413:	8b 00                	mov    (%rax),%eax
   11415:	39 85 d8 fb ff ff    	cmp    %eax,-0x428(%rbp)
   1141b:	7d 07                	jge    11424 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x757>
   1141d:	b8 01 00 00 00       	mov    $0x1,%eax
   11422:	eb 05                	jmp    11429 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x75c>
   11424:	b8 00 00 00 00       	mov    $0x0,%eax
   11429:	84 c0                	test   %al,%al
   1142b:	0f 85 49 fc ff ff    	jne    1107a <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x3ad>
   11431:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11438:	48 89 c7             	mov    %rax,%rdi
   1143b:	e8 50 8a ff ff       	call   9e90 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEE5closeEv@plt>
   11440:	48 8b 85 c0 fb ff ff 	mov    -0x440(%rbp),%rax
   11447:	48 8b 95 00 fc ff ff 	mov    -0x400(%rbp),%rdx
   1144e:	48 89 10             	mov    %rdx,(%rax)
   11451:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   11458:	8b 00                	mov    (%rax),%eax
   1145a:	39 85 d8 fb ff ff    	cmp    %eax,-0x428(%rbp)
   11460:	74 46                	je     114a8 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x7db>
   11462:	48 8b 85 b8 fb ff ff 	mov    -0x448(%rbp),%rax
   11469:	8b 08                	mov    (%rax),%ecx
   1146b:	48 8b 05 ae cd 00 00 	mov    0xcdae(%rip),%rax        # 1e220 <stderr@GLIBC_2.2.5>
   11472:	8b 95 d8 fb ff ff    	mov    -0x428(%rbp),%edx
   11478:	41 89 c9             	mov    %ecx,%r9d
   1147b:	41 89 d0             	mov    %edx,%r8d
   1147e:	b9 8c 00 00 00       	mov    $0x8c,%ecx
   11483:	48 8d 15 36 53 00 00 	lea    0x5336(%rip),%rdx        # 167c0 <_ZStL19piecewise_construct+0x40>
   1148a:	48 8d 35 67 54 00 00 	lea    0x5467(%rip),%rsi        # 168f8 <_ZStL19piecewise_construct+0x178>
   11491:	48 89 c7             	mov    %rax,%rdi
   11494:	b8 00 00 00 00       	mov    $0x0,%eax
   11499:	e8 b2 92 ff ff       	call   a750 <fprintf@plt>
   1149e:	bf 01 00 00 00       	mov    $0x1,%edi
   114a3:	e8 88 90 ff ff       	call   a530 <exit@plt>
   114a8:	48 8d 85 50 fc ff ff 	lea    -0x3b0(%rbp),%rax
   114af:	48 89 c7             	mov    %rax,%rdi
   114b2:	e8 39 8d ff ff       	call   a1f0 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEED1Ev@plt>
   114b7:	48 8d 85 30 fc ff ff 	lea    -0x3d0(%rbp),%rax
   114be:	48 89 c7             	mov    %rax,%rdi
   114c1:	e8 3a 8c ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   114c6:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   114cd:	48 89 c7             	mov    %rax,%rdi
   114d0:	e8 2b 8a ff ff       	call   9f00 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEED1Ev@plt>
   114d5:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   114d9:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   114e0:	00 00 
   114e2:	74 70                	je     11554 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x887>
   114e4:	eb 69                	jmp    1154f <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x882>
   114e6:	f3 0f 1e fa          	endbr64 
   114ea:	48 89 c3             	mov    %rax,%rbx
   114ed:	48 8d 85 50 fc ff ff 	lea    -0x3b0(%rbp),%rax
   114f4:	48 89 c7             	mov    %rax,%rdi
   114f7:	e8 04 8c ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   114fc:	eb 37                	jmp    11535 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x868>
   114fe:	f3 0f 1e fa          	endbr64 
   11502:	48 89 c3             	mov    %rax,%rbx
   11505:	48 8d 85 50 fc ff ff 	lea    -0x3b0(%rbp),%rax
   1150c:	48 89 c7             	mov    %rax,%rdi
   1150f:	e8 dc 8c ff ff       	call   a1f0 <_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEED1Ev@plt>
   11514:	eb 07                	jmp    1151d <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x850>
   11516:	f3 0f 1e fa          	endbr64 
   1151a:	48 89 c3             	mov    %rax,%rbx
   1151d:	48 8d 85 30 fc ff ff 	lea    -0x3d0(%rbp),%rax
   11524:	48 89 c7             	mov    %rax,%rdi
   11527:	e8 d4 8b ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   1152c:	eb 07                	jmp    11535 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi+0x868>
   1152e:	f3 0f 1e fa          	endbr64 
   11532:	48 89 c3             	mov    %rax,%rbx
   11535:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   1153c:	48 89 c7             	mov    %rax,%rdi
   1153f:	e8 bc 89 ff ff       	call   9f00 <_ZNSt14basic_ifstreamIcSt11char_traitsIcEED1Ev@plt>
   11544:	48 89 d8             	mov    %rbx,%rax
   11547:	48 89 c7             	mov    %rax,%rdi
   1154a:	e8 31 94 ff ff       	call   a980 <_Unwind_Resume@plt>
   1154f:	e8 7c 8f ff ff       	call   a4d0 <__stack_chk_fail@plt>
   11554:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   11558:	c9                   	leave  
   11559:	c3                   	ret    

000000000001155a <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti>:
   1155a:	f3 0f 1e fa          	endbr64 
   1155e:	55                   	push   %rbp
   1155f:	48 89 e5             	mov    %rsp,%rbp
   11562:	53                   	push   %rbx
   11563:	48 81 ec 78 02 00 00 	sub    $0x278,%rsp
   1156a:	48 89 bd 98 fd ff ff 	mov    %rdi,-0x268(%rbp)
   11571:	48 89 b5 90 fd ff ff 	mov    %rsi,-0x270(%rbp)
   11578:	89 95 8c fd ff ff    	mov    %edx,-0x274(%rbp)
   1157e:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   11585:	00 00 
   11587:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   1158b:	31 c0                	xor    %eax,%eax
   1158d:	be 20 00 00 00       	mov    $0x20,%esi
   11592:	bf 10 00 00 00       	mov    $0x10,%edi
   11597:	e8 a3 0c 00 00       	call   1223f <_ZStorSt13_Ios_OpenmodeS_>
   1159c:	89 c2                	mov    %eax,%edx
   1159e:	48 8b 8d 98 fd ff ff 	mov    -0x268(%rbp),%rcx
   115a5:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   115ac:	48 89 ce             	mov    %rcx,%rsi
   115af:	48 89 c7             	mov    %rax,%rdi
   115b2:	e8 f9 94 ff ff       	call   aab0 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEEC1ERKNSt7__cxx1112basic_stringIcS1_SaIcEEESt13_Ios_Openmode@plt>
   115b7:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   115be:	48 89 c7             	mov    %rax,%rdi
   115c1:	e8 6a 8c ff ff       	call   a230 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEE7is_openEv@plt>
   115c6:	83 f0 01             	xor    $0x1,%eax
   115c9:	84 c0                	test   %al,%al
   115cb:	74 63                	je     11630 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0xd6>
   115cd:	48 8b 95 98 fd ff ff 	mov    -0x268(%rbp),%rdx
   115d4:	48 8d 85 c0 fd ff ff 	lea    -0x240(%rbp),%rax
   115db:	48 89 d6             	mov    %rdx,%rsi
   115de:	48 89 c7             	mov    %rax,%rdi
   115e1:	e8 da 89 ff ff       	call   9fc0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
   115e6:	48 8b 05 33 cc 00 00 	mov    0xcc33(%rip),%rax        # 1e220 <stderr@GLIBC_2.2.5>
   115ed:	48 8d 95 c0 fd ff ff 	lea    -0x240(%rbp),%rdx
   115f4:	49 89 d0             	mov    %rdx,%r8
   115f7:	b9 94 00 00 00       	mov    $0x94,%ecx
   115fc:	48 8d 15 bd 51 00 00 	lea    0x51bd(%rip),%rdx        # 167c0 <_ZStL19piecewise_construct+0x40>
   11603:	48 8d 35 46 52 00 00 	lea    0x5246(%rip),%rsi        # 16850 <_ZStL19piecewise_construct+0xd0>
   1160a:	48 89 c7             	mov    %rax,%rdi
   1160d:	b8 00 00 00 00       	mov    $0x0,%eax
   11612:	e8 39 91 ff ff       	call   a750 <fprintf@plt>
   11617:	48 8d 85 c0 fd ff ff 	lea    -0x240(%rbp),%rax
   1161e:	48 89 c7             	mov    %rax,%rdi
   11621:	e8 da 8a ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   11626:	bf 01 00 00 00       	mov    $0x1,%edi
   1162b:	e8 00 8f ff ff       	call   a530 <exit@plt>
   11630:	48 8b 85 90 fd ff ff 	mov    -0x270(%rbp),%rax
   11637:	8b 40 30             	mov    0x30(%rax),%eax
   1163a:	89 85 b4 fd ff ff    	mov    %eax,-0x24c(%rbp)
   11640:	8b 95 8c fd ff ff    	mov    -0x274(%rbp),%edx
   11646:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   1164d:	89 d6                	mov    %edx,%esi
   1164f:	48 89 c7             	mov    %rax,%rdi
   11652:	e8 f9 92 ff ff       	call   a950 <_ZNSolsEi@plt>
   11657:	48 89 c2             	mov    %rax,%rdx
   1165a:	48 8d 05 81 52 00 00 	lea    0x5281(%rip),%rax        # 168e2 <_ZStL19piecewise_construct+0x162>
   11661:	48 89 c6             	mov    %rax,%rsi
   11664:	48 89 d7             	mov    %rdx,%rdi
   11667:	e8 54 8d ff ff       	call   a3c0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
   1166c:	48 89 c2             	mov    %rax,%rdx
   1166f:	8b 85 b4 fd ff ff    	mov    -0x24c(%rbp),%eax
   11675:	89 c6                	mov    %eax,%esi
   11677:	48 89 d7             	mov    %rdx,%rdi
   1167a:	e8 d1 92 ff ff       	call   a950 <_ZNSolsEi@plt>
   1167f:	48 8b 15 4a b9 00 00 	mov    0xb94a(%rip),%rdx        # 1cfd0 <_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_@GLIBCXX_3.4>
   11686:	48 89 d6             	mov    %rdx,%rsi
   11689:	48 89 c7             	mov    %rax,%rdi
   1168c:	e8 6f 8d ff ff       	call   a400 <_ZNSolsEPFRSoS_E@plt>
   11691:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11698:	48 89 c7             	mov    %rax,%rdi
   1169b:	e8 d0 8a ff ff       	call   a170 <_ZNSo5flushEv@plt>
   116a0:	48 8b 85 90 fd ff ff 	mov    -0x270(%rbp),%rax
   116a7:	48 89 85 b8 fd ff ff 	mov    %rax,-0x248(%rbp)
   116ae:	c7 85 a0 fd ff ff 00 	movl   $0x0,-0x260(%rbp)
   116b5:	00 00 00 
   116b8:	e9 da 02 00 00       	jmp    11997 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x43d>
   116bd:	8b 85 a0 fd ff ff    	mov    -0x260(%rbp),%eax
   116c3:	48 98                	cltq   
   116c5:	48 c1 e0 06          	shl    $0x6,%rax
   116c9:	48 89 c2             	mov    %rax,%rdx
   116cc:	48 8b 85 90 fd ff ff 	mov    -0x270(%rbp),%rax
   116d3:	48 01 d0             	add    %rdx,%rax
   116d6:	48 89 85 b8 fd ff ff 	mov    %rax,-0x248(%rbp)
   116dd:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   116e4:	48 8d 15 d1 0b 00 00 	lea    0xbd1(%rip),%rdx        # 122bc <_ZSt10scientificRSt8ios_base>
   116eb:	48 89 d6             	mov    %rdx,%rsi
   116ee:	48 89 c7             	mov    %rax,%rdi
   116f1:	e8 7a 88 ff ff       	call   9f70 <_ZNSolsEPFRSt8ios_baseS0_E@plt>
   116f6:	48 89 c2             	mov    %rax,%rdx
   116f9:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11700:	48 8b 00             	mov    (%rax),%rax
   11703:	48 89 c6             	mov    %rax,%rsi
   11706:	48 89 d7             	mov    %rdx,%rdi
   11709:	e8 c2 87 ff ff       	call   9ed0 <_ZNSolsEm@plt>
   1170e:	48 89 c2             	mov    %rax,%rdx
   11711:	48 8d 05 ca 51 00 00 	lea    0x51ca(%rip),%rax        # 168e2 <_ZStL19piecewise_construct+0x162>
   11718:	48 89 c6             	mov    %rax,%rsi
   1171b:	48 89 d7             	mov    %rdx,%rdi
   1171e:	e8 9d 8c ff ff       	call   a3c0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
   11723:	c7 85 a4 fd ff ff 00 	movl   $0x0,-0x25c(%rbp)
   1172a:	00 00 00 
   1172d:	eb 62                	jmp    11791 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x237>
   1172f:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11736:	48 8d 15 7f 0b 00 00 	lea    0xb7f(%rip),%rdx        # 122bc <_ZSt10scientificRSt8ios_base>
   1173d:	48 89 d6             	mov    %rdx,%rsi
   11740:	48 89 c7             	mov    %rax,%rdi
   11743:	e8 28 88 ff ff       	call   9f70 <_ZNSolsEPFRSt8ios_baseS0_E@plt>
   11748:	48 89 c2             	mov    %rax,%rdx
   1174b:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11752:	48 8b 48 08          	mov    0x8(%rax),%rcx
   11756:	8b 85 a4 fd ff ff    	mov    -0x25c(%rbp),%eax
   1175c:	48 98                	cltq   
   1175e:	48 c1 e0 03          	shl    $0x3,%rax
   11762:	48 01 c8             	add    %rcx,%rax
   11765:	48 8b 00             	mov    (%rax),%rax
   11768:	66 48 0f 6e c0       	movq   %rax,%xmm0
   1176d:	48 89 d7             	mov    %rdx,%rdi
   11770:	e8 cb 92 ff ff       	call   aa40 <_ZNSolsEd@plt>
   11775:	48 89 c2             	mov    %rax,%rdx
   11778:	48 8d 05 63 51 00 00 	lea    0x5163(%rip),%rax        # 168e2 <_ZStL19piecewise_construct+0x162>
   1177f:	48 89 c6             	mov    %rax,%rsi
   11782:	48 89 d7             	mov    %rdx,%rdi
   11785:	e8 36 8c ff ff       	call   a3c0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
   1178a:	83 85 a4 fd ff ff 01 	addl   $0x1,-0x25c(%rbp)
   11791:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11798:	8b 40 30             	mov    0x30(%rax),%eax
   1179b:	39 85 a4 fd ff ff    	cmp    %eax,-0x25c(%rbp)
   117a1:	7c 8c                	jl     1172f <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x1d5>
   117a3:	c7 85 a8 fd ff ff 00 	movl   $0x0,-0x258(%rbp)
   117aa:	00 00 00 
   117ad:	eb 62                	jmp    11811 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x2b7>
   117af:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   117b6:	48 8d 15 ff 0a 00 00 	lea    0xaff(%rip),%rdx        # 122bc <_ZSt10scientificRSt8ios_base>
   117bd:	48 89 d6             	mov    %rdx,%rsi
   117c0:	48 89 c7             	mov    %rax,%rdi
   117c3:	e8 a8 87 ff ff       	call   9f70 <_ZNSolsEPFRSt8ios_baseS0_E@plt>
   117c8:	48 89 c2             	mov    %rax,%rdx
   117cb:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   117d2:	48 8b 48 10          	mov    0x10(%rax),%rcx
   117d6:	8b 85 a8 fd ff ff    	mov    -0x258(%rbp),%eax
   117dc:	48 98                	cltq   
   117de:	48 c1 e0 03          	shl    $0x3,%rax
   117e2:	48 01 c8             	add    %rcx,%rax
   117e5:	48 8b 00             	mov    (%rax),%rax
   117e8:	66 48 0f 6e c0       	movq   %rax,%xmm0
   117ed:	48 89 d7             	mov    %rdx,%rdi
   117f0:	e8 4b 92 ff ff       	call   aa40 <_ZNSolsEd@plt>
   117f5:	48 89 c2             	mov    %rax,%rdx
   117f8:	48 8d 05 e3 50 00 00 	lea    0x50e3(%rip),%rax        # 168e2 <_ZStL19piecewise_construct+0x162>
   117ff:	48 89 c6             	mov    %rax,%rsi
   11802:	48 89 d7             	mov    %rdx,%rdi
   11805:	e8 b6 8b ff ff       	call   a3c0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
   1180a:	83 85 a8 fd ff ff 01 	addl   $0x1,-0x258(%rbp)
   11811:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11818:	8b 40 30             	mov    0x30(%rax),%eax
   1181b:	39 85 a8 fd ff ff    	cmp    %eax,-0x258(%rbp)
   11821:	7c 8c                	jl     117af <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x255>
   11823:	c7 85 ac fd ff ff 00 	movl   $0x0,-0x254(%rbp)
   1182a:	00 00 00 
   1182d:	eb 62                	jmp    11891 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x337>
   1182f:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11836:	48 8d 15 7f 0a 00 00 	lea    0xa7f(%rip),%rdx        # 122bc <_ZSt10scientificRSt8ios_base>
   1183d:	48 89 d6             	mov    %rdx,%rsi
   11840:	48 89 c7             	mov    %rax,%rdi
   11843:	e8 28 87 ff ff       	call   9f70 <_ZNSolsEPFRSt8ios_baseS0_E@plt>
   11848:	48 89 c2             	mov    %rax,%rdx
   1184b:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11852:	48 8b 48 18          	mov    0x18(%rax),%rcx
   11856:	8b 85 ac fd ff ff    	mov    -0x254(%rbp),%eax
   1185c:	48 98                	cltq   
   1185e:	48 c1 e0 03          	shl    $0x3,%rax
   11862:	48 01 c8             	add    %rcx,%rax
   11865:	48 8b 00             	mov    (%rax),%rax
   11868:	66 48 0f 6e c0       	movq   %rax,%xmm0
   1186d:	48 89 d7             	mov    %rdx,%rdi
   11870:	e8 cb 91 ff ff       	call   aa40 <_ZNSolsEd@plt>
   11875:	48 89 c2             	mov    %rax,%rdx
   11878:	48 8d 05 63 50 00 00 	lea    0x5063(%rip),%rax        # 168e2 <_ZStL19piecewise_construct+0x162>
   1187f:	48 89 c6             	mov    %rax,%rsi
   11882:	48 89 d7             	mov    %rdx,%rdi
   11885:	e8 36 8b ff ff       	call   a3c0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
   1188a:	83 85 ac fd ff ff 01 	addl   $0x1,-0x254(%rbp)
   11891:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11898:	8b 40 30             	mov    0x30(%rax),%eax
   1189b:	39 85 ac fd ff ff    	cmp    %eax,-0x254(%rbp)
   118a1:	7c 8c                	jl     1182f <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x2d5>
   118a3:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   118aa:	8b 50 34             	mov    0x34(%rax),%edx
   118ad:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   118b4:	89 d6                	mov    %edx,%esi
   118b6:	48 89 c7             	mov    %rax,%rdi
   118b9:	e8 92 90 ff ff       	call   a950 <_ZNSolsEi@plt>
   118be:	48 89 c2             	mov    %rax,%rdx
   118c1:	48 8d 05 1a 50 00 00 	lea    0x501a(%rip),%rax        # 168e2 <_ZStL19piecewise_construct+0x162>
   118c8:	48 89 c6             	mov    %rax,%rsi
   118cb:	48 89 d7             	mov    %rdx,%rdi
   118ce:	e8 ed 8a ff ff       	call   a3c0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
   118d3:	c7 85 b0 fd ff ff 00 	movl   $0x0,-0x250(%rbp)
   118da:	00 00 00 
   118dd:	e9 98 00 00 00       	jmp    1197a <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x420>
   118e2:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   118e9:	48 8d 15 cc 09 00 00 	lea    0x9cc(%rip),%rdx        # 122bc <_ZSt10scientificRSt8ios_base>
   118f0:	48 89 d6             	mov    %rdx,%rsi
   118f3:	48 89 c7             	mov    %rax,%rdi
   118f6:	e8 75 86 ff ff       	call   9f70 <_ZNSolsEPFRSt8ios_baseS0_E@plt>
   118fb:	48 89 c2             	mov    %rax,%rdx
   118fe:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11905:	48 8b 48 20          	mov    0x20(%rax),%rcx
   11909:	8b 85 b0 fd ff ff    	mov    -0x250(%rbp),%eax
   1190f:	48 98                	cltq   
   11911:	48 c1 e0 03          	shl    $0x3,%rax
   11915:	48 01 c8             	add    %rcx,%rax
   11918:	48 8b 00             	mov    (%rax),%rax
   1191b:	66 48 0f 6e c0       	movq   %rax,%xmm0
   11920:	48 89 d7             	mov    %rdx,%rdi
   11923:	e8 18 91 ff ff       	call   aa40 <_ZNSolsEd@plt>
   11928:	8b 85 b0 fd ff ff    	mov    -0x250(%rbp),%eax
   1192e:	8d 50 01             	lea    0x1(%rax),%edx
   11931:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11938:	8b 40 34             	mov    0x34(%rax),%eax
   1193b:	39 c2                	cmp    %eax,%edx
   1193d:	74 1b                	je     1195a <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x400>
   1193f:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11946:	48 8d 15 95 4f 00 00 	lea    0x4f95(%rip),%rdx        # 168e2 <_ZStL19piecewise_construct+0x162>
   1194d:	48 89 d6             	mov    %rdx,%rsi
   11950:	48 89 c7             	mov    %rax,%rdi
   11953:	e8 68 8a ff ff       	call   a3c0 <_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc@plt>
   11958:	eb 19                	jmp    11973 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x419>
   1195a:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   11961:	48 8b 15 68 b6 00 00 	mov    0xb668(%rip),%rdx        # 1cfd0 <_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_@GLIBCXX_3.4>
   11968:	48 89 d6             	mov    %rdx,%rsi
   1196b:	48 89 c7             	mov    %rax,%rdi
   1196e:	e8 8d 8a ff ff       	call   a400 <_ZNSolsEPFRSoS_E@plt>
   11973:	83 85 b0 fd ff ff 01 	addl   $0x1,-0x250(%rbp)
   1197a:	48 8b 85 b8 fd ff ff 	mov    -0x248(%rbp),%rax
   11981:	8b 40 34             	mov    0x34(%rax),%eax
   11984:	39 85 b0 fd ff ff    	cmp    %eax,-0x250(%rbp)
   1198a:	0f 8c 52 ff ff ff    	jl     118e2 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x388>
   11990:	83 85 a0 fd ff ff 01 	addl   $0x1,-0x260(%rbp)
   11997:	8b 85 a0 fd ff ff    	mov    -0x260(%rbp),%eax
   1199d:	3b 85 8c fd ff ff    	cmp    -0x274(%rbp),%eax
   119a3:	0f 8c 14 fd ff ff    	jl     116bd <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x163>
   119a9:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   119b0:	48 89 c7             	mov    %rax,%rdi
   119b3:	e8 58 8d ff ff       	call   a710 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEE5closeEv@plt>
   119b8:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   119bf:	48 89 c7             	mov    %rax,%rdi
   119c2:	e8 a9 8d ff ff       	call   a770 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEED1Ev@plt>
   119c7:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   119cb:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   119d2:	00 00 
   119d4:	74 40                	je     11a16 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x4bc>
   119d6:	eb 39                	jmp    11a11 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x4b7>
   119d8:	f3 0f 1e fa          	endbr64 
   119dc:	48 89 c3             	mov    %rax,%rbx
   119df:	48 8d 85 c0 fd ff ff 	lea    -0x240(%rbp),%rax
   119e6:	48 89 c7             	mov    %rax,%rdi
   119e9:	e8 12 87 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   119ee:	eb 07                	jmp    119f7 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti+0x49d>
   119f0:	f3 0f 1e fa          	endbr64 
   119f4:	48 89 c3             	mov    %rax,%rbx
   119f7:	48 8d 85 e0 fd ff ff 	lea    -0x220(%rbp),%rax
   119fe:	48 89 c7             	mov    %rax,%rdi
   11a01:	e8 6a 8d ff ff       	call   a770 <_ZNSt14basic_ofstreamIcSt11char_traitsIcEED1Ev@plt>
   11a06:	48 89 d8             	mov    %rbx,%rax
   11a09:	48 89 c7             	mov    %rax,%rdi
   11a0c:	e8 6f 8f ff ff       	call   a980 <_Unwind_Resume@plt>
   11a11:	e8 ba 8a ff ff       	call   a4d0 <__stack_chk_fail@plt>
   11a16:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   11a1a:	c9                   	leave  
   11a1b:	c3                   	ret    

0000000000011a1c <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t>:
   11a1c:	f3 0f 1e fa          	endbr64 
   11a20:	55                   	push   %rbp
   11a21:	48 89 e5             	mov    %rsp,%rbp
   11a24:	53                   	push   %rbx
   11a25:	48 83 ec 58          	sub    $0x58,%rsp
   11a29:	48 89 7d a8          	mov    %rdi,-0x58(%rbp)
   11a2d:	48 89 75 a0          	mov    %rsi,-0x60(%rbp)
   11a31:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   11a38:	00 00 
   11a3a:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   11a3e:	31 c0                	xor    %eax,%eax
   11a40:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   11a44:	48 8b 00             	mov    (%rax),%rax
   11a47:	48 85 c0             	test   %rax,%rax
   11a4a:	75 16                	jne    11a62 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t+0x46>
   11a4c:	bf 10 00 00 00       	mov    $0x10,%edi
   11a51:	e8 ca 8b ff ff       	call   a620 <malloc@plt>
   11a56:	48 89 c2             	mov    %rax,%rdx
   11a59:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   11a5d:	48 89 10             	mov    %rdx,(%rax)
   11a60:	eb 36                	jmp    11a98 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t+0x7c>
   11a62:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   11a66:	48 8b 00             	mov    (%rax),%rax
   11a69:	8b 40 08             	mov    0x8(%rax),%eax
   11a6c:	85 c0                	test   %eax,%eax
   11a6e:	74 28                	je     11a98 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t+0x7c>
   11a70:	48 8d 05 b1 4e 00 00 	lea    0x4eb1(%rip),%rax        # 16928 <_ZStL19piecewise_construct+0x1a8>
   11a77:	48 89 c1             	mov    %rax,%rcx
   11a7a:	ba b9 00 00 00       	mov    $0xb9,%edx
   11a7f:	48 8d 05 3a 4d 00 00 	lea    0x4d3a(%rip),%rax        # 167c0 <_ZStL19piecewise_construct+0x40>
   11a86:	48 89 c6             	mov    %rax,%rsi
   11a89:	48 8d 05 d0 4e 00 00 	lea    0x4ed0(%rip),%rax        # 16960 <_ZStL19piecewise_construct+0x1e0>
   11a90:	48 89 c7             	mov    %rax,%rdi
   11a93:	e8 98 86 ff ff       	call   a130 <__assert_fail@plt>
   11a98:	48 8b 55 a8          	mov    -0x58(%rbp),%rdx
   11a9c:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11aa0:	48 89 d6             	mov    %rdx,%rsi
   11aa3:	48 89 c7             	mov    %rax,%rdi
   11aa6:	e8 15 85 ff ff       	call   9fc0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
   11aab:	48 8d 55 b4          	lea    -0x4c(%rbp),%rdx
   11aaf:	48 8d 4d b8          	lea    -0x48(%rbp),%rcx
   11ab3:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11ab7:	48 89 ce             	mov    %rcx,%rsi
   11aba:	48 89 c7             	mov    %rax,%rdi
   11abd:	e8 0b f2 ff ff       	call   10ccd <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP10particle_tPi>
   11ac2:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11ac6:	48 89 c7             	mov    %rax,%rdi
   11ac9:	e8 32 86 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   11ace:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   11ad2:	48 8b 00             	mov    (%rax),%rax
   11ad5:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
   11ad9:	48 89 10             	mov    %rdx,(%rax)
   11adc:	48 8b 45 a0          	mov    -0x60(%rbp),%rax
   11ae0:	48 8b 00             	mov    (%rax),%rax
   11ae3:	8b 55 b4             	mov    -0x4c(%rbp),%edx
   11ae6:	89 50 08             	mov    %edx,0x8(%rax)
   11ae9:	90                   	nop
   11aea:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   11aee:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   11af5:	00 00 
   11af7:	74 25                	je     11b1e <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t+0x102>
   11af9:	eb 1e                	jmp    11b19 <_Z14load_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEPP17chunk_particles_t+0xfd>
   11afb:	f3 0f 1e fa          	endbr64 
   11aff:	48 89 c3             	mov    %rax,%rbx
   11b02:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11b06:	48 89 c7             	mov    %rax,%rdi
   11b09:	e8 f2 85 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   11b0e:	48 89 d8             	mov    %rbx,%rax
   11b11:	48 89 c7             	mov    %rax,%rdi
   11b14:	e8 67 8e ff ff       	call   a980 <_Unwind_Resume@plt>
   11b19:	e8 b2 89 ff ff       	call   a4d0 <__stack_chk_fail@plt>
   11b1e:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   11b22:	c9                   	leave  
   11b23:	c3                   	ret    

0000000000011b24 <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t>:
   11b24:	f3 0f 1e fa          	endbr64 
   11b28:	55                   	push   %rbp
   11b29:	48 89 e5             	mov    %rsp,%rbp
   11b2c:	41 54                	push   %r12
   11b2e:	53                   	push   %rbx
   11b2f:	48 83 ec 40          	sub    $0x40,%rsp
   11b33:	48 89 7d b8          	mov    %rdi,-0x48(%rbp)
   11b37:	48 89 75 b0          	mov    %rsi,-0x50(%rbp)
   11b3b:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   11b42:	00 00 
   11b44:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   11b48:	31 c0                	xor    %eax,%eax
   11b4a:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   11b4e:	44 8b 60 08          	mov    0x8(%rax),%r12d
   11b52:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   11b56:	48 8b 18             	mov    (%rax),%rbx
   11b59:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
   11b5d:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11b61:	48 89 d6             	mov    %rdx,%rsi
   11b64:	48 89 c7             	mov    %rax,%rdi
   11b67:	e8 54 84 ff ff       	call   9fc0 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1ERKS4_@plt>
   11b6c:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11b70:	44 89 e2             	mov    %r12d,%edx
   11b73:	48 89 de             	mov    %rbx,%rsi
   11b76:	48 89 c7             	mov    %rax,%rdi
   11b79:	e8 dc f9 ff ff       	call   1155a <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP10particle_ti>
   11b7e:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11b82:	48 89 c7             	mov    %rax,%rdi
   11b85:	e8 76 85 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   11b8a:	90                   	nop
   11b8b:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   11b8f:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   11b96:	00 00 
   11b98:	74 25                	je     11bbf <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t+0x9b>
   11b9a:	eb 1e                	jmp    11bba <_Z15store_particlesNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEP17chunk_particles_t+0x96>
   11b9c:	f3 0f 1e fa          	endbr64 
   11ba0:	48 89 c3             	mov    %rax,%rbx
   11ba3:	48 8d 45 c0          	lea    -0x40(%rbp),%rax
   11ba7:	48 89 c7             	mov    %rax,%rdi
   11baa:	e8 51 85 ff ff       	call   a100 <_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEED1Ev@plt>
   11baf:	48 89 d8             	mov    %rbx,%rax
   11bb2:	48 89 c7             	mov    %rax,%rdi
   11bb5:	e8 c6 8d ff ff       	call   a980 <_Unwind_Resume@plt>
   11bba:	e8 11 89 ff ff       	call   a4d0 <__stack_chk_fail@plt>
   11bbf:	48 83 c4 40          	add    $0x40,%rsp
   11bc3:	5b                   	pop    %rbx
   11bc4:	41 5c                	pop    %r12
   11bc6:	5d                   	pop    %rbp
   11bc7:	c3                   	ret    

0000000000011bc8 <_Z14free_particlesP10particle_t>:
   11bc8:	f3 0f 1e fa          	endbr64 
   11bcc:	55                   	push   %rbp
   11bcd:	48 89 e5             	mov    %rsp,%rbp
   11bd0:	48 83 ec 10          	sub    $0x10,%rsp
   11bd4:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   11bd8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11bdc:	48 8b 40 08          	mov    0x8(%rax),%rax
   11be0:	48 89 c7             	mov    %rax,%rdi
   11be3:	e8 08 89 ff ff       	call   a4f0 <free@plt>
   11be8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11bec:	48 8b 40 10          	mov    0x10(%rax),%rax
   11bf0:	48 89 c7             	mov    %rax,%rdi
   11bf3:	e8 f8 88 ff ff       	call   a4f0 <free@plt>
   11bf8:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11bfc:	48 8b 40 18          	mov    0x18(%rax),%rax
   11c00:	48 89 c7             	mov    %rax,%rdi
   11c03:	e8 e8 88 ff ff       	call   a4f0 <free@plt>
   11c08:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11c0c:	48 8b 40 20          	mov    0x20(%rax),%rax
   11c10:	48 89 c7             	mov    %rax,%rdi
   11c13:	e8 d8 88 ff ff       	call   a4f0 <free@plt>
   11c18:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11c1c:	48 89 c7             	mov    %rax,%rdi
   11c1f:	e8 cc 88 ff ff       	call   a4f0 <free@plt>
   11c24:	90                   	nop
   11c25:	c9                   	leave  
   11c26:	c3                   	ret    

0000000000011c27 <_Z14free_particlesP17chunk_particles_t>:
   11c27:	f3 0f 1e fa          	endbr64 
   11c2b:	55                   	push   %rbp
   11c2c:	48 89 e5             	mov    %rsp,%rbp
   11c2f:	48 83 ec 10          	sub    $0x10,%rsp
   11c33:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   11c37:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11c3b:	48 8b 00             	mov    (%rax),%rax
   11c3e:	48 89 c7             	mov    %rax,%rdi
   11c41:	e8 82 ff ff ff       	call   11bc8 <_Z14free_particlesP10particle_t>
   11c46:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11c4a:	48 c7 00 00 00 00 00 	movq   $0x0,(%rax)
   11c51:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11c55:	c7 40 08 00 00 00 00 	movl   $0x0,0x8(%rax)
   11c5c:	90                   	nop
   11c5d:	c9                   	leave  
   11c5e:	c3                   	ret    

0000000000011c5f <_Z18clearAccelerationsP10particle_ti>:
   11c5f:	f3 0f 1e fa          	endbr64 
   11c63:	55                   	push   %rbp
   11c64:	48 89 e5             	mov    %rsp,%rbp
   11c67:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   11c6b:	89 75 d4             	mov    %esi,-0x2c(%rbp)
   11c6e:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   11c72:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   11c76:	c7 45 e8 00 00 00 00 	movl   $0x0,-0x18(%rbp)
   11c7d:	eb 42                	jmp    11cc1 <_Z18clearAccelerationsP10particle_ti+0x62>
   11c7f:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   11c83:	48 8b 40 18          	mov    0x18(%rax),%rax
   11c87:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   11c8b:	c7 45 ec 00 00 00 00 	movl   $0x0,-0x14(%rbp)
   11c92:	eb 18                	jmp    11cac <_Z18clearAccelerationsP10particle_ti+0x4d>
   11c94:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11c98:	48 8d 50 08          	lea    0x8(%rax),%rdx
   11c9c:	48 89 55 f8          	mov    %rdx,-0x8(%rbp)
   11ca0:	66 0f ef c0          	pxor   %xmm0,%xmm0
   11ca4:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   11ca8:	83 45 ec 01          	addl   $0x1,-0x14(%rbp)
   11cac:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   11cb0:	8b 40 30             	mov    0x30(%rax),%eax
   11cb3:	39 45 ec             	cmp    %eax,-0x14(%rbp)
   11cb6:	7c dc                	jl     11c94 <_Z18clearAccelerationsP10particle_ti+0x35>
   11cb8:	48 83 45 f0 40       	addq   $0x40,-0x10(%rbp)
   11cbd:	83 45 e8 01          	addl   $0x1,-0x18(%rbp)
   11cc1:	8b 45 e8             	mov    -0x18(%rbp),%eax
   11cc4:	3b 45 d4             	cmp    -0x2c(%rbp),%eax
   11cc7:	7c b6                	jl     11c7f <_Z18clearAccelerationsP10particle_ti+0x20>
   11cc9:	90                   	nop
   11cca:	90                   	nop
   11ccb:	5d                   	pop    %rbp
   11ccc:	c3                   	ret    

0000000000011ccd <_Z18clearAccelerationsP17chunk_particles_t>:
   11ccd:	f3 0f 1e fa          	endbr64 
   11cd1:	55                   	push   %rbp
   11cd2:	48 89 e5             	mov    %rsp,%rbp
   11cd5:	48 83 ec 08          	sub    $0x8,%rsp
   11cd9:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   11cdd:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11ce1:	8b 50 08             	mov    0x8(%rax),%edx
   11ce4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11ce8:	48 8b 00             	mov    (%rax),%rax
   11ceb:	89 d6                	mov    %edx,%esi
   11ced:	48 89 c7             	mov    %rax,%rdi
   11cf0:	e8 6a ff ff ff       	call   11c5f <_Z18clearAccelerationsP10particle_ti>
   11cf5:	90                   	nop
   11cf6:	c9                   	leave  
   11cf7:	c3                   	ret    

0000000000011cf8 <_Z16update_particlesP10particle_tid>:
   11cf8:	f3 0f 1e fa          	endbr64 
   11cfc:	55                   	push   %rbp
   11cfd:	48 89 e5             	mov    %rsp,%rbp
   11d00:	48 83 ec 30          	sub    $0x30,%rsp
   11d04:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   11d08:	89 75 e4             	mov    %esi,-0x1c(%rbp)
   11d0b:	f2 0f 11 45 d8       	movsd  %xmm0,-0x28(%rbp)
   11d10:	8b 55 e4             	mov    -0x1c(%rbp),%edx
   11d13:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   11d17:	89 d6                	mov    %edx,%esi
   11d19:	48 89 c7             	mov    %rax,%rdi
   11d1c:	e8 3e ff ff ff       	call   11c5f <_Z18clearAccelerationsP10particle_ti>
   11d21:	c7 45 f4 00 00 00 00 	movl   $0x0,-0xc(%rbp)
   11d28:	eb 4a                	jmp    11d74 <_Z16update_particlesP10particle_tid+0x7c>
   11d2a:	c7 45 f8 00 00 00 00 	movl   $0x0,-0x8(%rbp)
   11d31:	eb 35                	jmp    11d68 <_Z16update_particlesP10particle_tid+0x70>
   11d33:	8b 45 f8             	mov    -0x8(%rbp),%eax
   11d36:	48 98                	cltq   
   11d38:	48 c1 e0 06          	shl    $0x6,%rax
   11d3c:	48 89 c2             	mov    %rax,%rdx
   11d3f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   11d43:	48 01 c2             	add    %rax,%rdx
   11d46:	8b 45 f4             	mov    -0xc(%rbp),%eax
   11d49:	48 98                	cltq   
   11d4b:	48 c1 e0 06          	shl    $0x6,%rax
   11d4f:	48 89 c1             	mov    %rax,%rcx
   11d52:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   11d56:	48 01 c8             	add    %rcx,%rax
   11d59:	48 89 d6             	mov    %rdx,%rsi
   11d5c:	48 89 c7             	mov    %rax,%rdi
   11d5f:	e8 bd ed ff ff       	call   10b21 <_Z18updateAccelerationP10particle_tS0_>
   11d64:	83 45 f8 01          	addl   $0x1,-0x8(%rbp)
   11d68:	8b 45 f8             	mov    -0x8(%rbp),%eax
   11d6b:	3b 45 f4             	cmp    -0xc(%rbp),%eax
   11d6e:	7c c3                	jl     11d33 <_Z16update_particlesP10particle_tid+0x3b>
   11d70:	83 45 f4 01          	addl   $0x1,-0xc(%rbp)
   11d74:	8b 45 f4             	mov    -0xc(%rbp),%eax
   11d77:	3b 45 e4             	cmp    -0x1c(%rbp),%eax
   11d7a:	7c ae                	jl     11d2a <_Z16update_particlesP10particle_tid+0x32>
   11d7c:	c7 45 fc 00 00 00 00 	movl   $0x0,-0x4(%rbp)
   11d83:	eb 4c                	jmp    11dd1 <_Z16update_particlesP10particle_tid+0xd9>
   11d85:	8b 45 fc             	mov    -0x4(%rbp),%eax
   11d88:	48 98                	cltq   
   11d8a:	48 c1 e0 06          	shl    $0x6,%rax
   11d8e:	48 89 c2             	mov    %rax,%rdx
   11d91:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   11d95:	48 01 c2             	add    %rax,%rdx
   11d98:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   11d9c:	66 48 0f 6e c0       	movq   %rax,%xmm0
   11da1:	48 89 d7             	mov    %rdx,%rdi
   11da4:	e8 6a ec ff ff       	call   10a13 <_Z14updateVelocityP10particle_td>
   11da9:	8b 45 fc             	mov    -0x4(%rbp),%eax
   11dac:	48 98                	cltq   
   11dae:	48 c1 e0 06          	shl    $0x6,%rax
   11db2:	48 89 c2             	mov    %rax,%rdx
   11db5:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   11db9:	48 01 c2             	add    %rax,%rdx
   11dbc:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   11dc0:	66 48 0f 6e c0       	movq   %rax,%xmm0
   11dc5:	48 89 d7             	mov    %rdx,%rdi
   11dc8:	e8 cd ec ff ff       	call   10a9a <_Z14updatePositionP10particle_td>
   11dcd:	83 45 fc 01          	addl   $0x1,-0x4(%rbp)
   11dd1:	8b 45 fc             	mov    -0x4(%rbp),%eax
   11dd4:	3b 45 e4             	cmp    -0x1c(%rbp),%eax
   11dd7:	7c ac                	jl     11d85 <_Z16update_particlesP10particle_tid+0x8d>
   11dd9:	90                   	nop
   11dda:	90                   	nop
   11ddb:	c9                   	leave  
   11ddc:	c3                   	ret    

0000000000011ddd <_Z16update_particlesP17chunk_particles_td>:
   11ddd:	f3 0f 1e fa          	endbr64 
   11de1:	55                   	push   %rbp
   11de2:	48 89 e5             	mov    %rsp,%rbp
   11de5:	48 83 ec 10          	sub    $0x10,%rsp
   11de9:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   11ded:	f2 0f 11 45 f0       	movsd  %xmm0,-0x10(%rbp)
   11df2:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11df6:	8b 50 08             	mov    0x8(%rax),%edx
   11df9:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   11dfd:	48 8b 00             	mov    (%rax),%rax
   11e00:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
   11e04:	66 48 0f 6e c1       	movq   %rcx,%xmm0
   11e09:	89 d6                	mov    %edx,%esi
   11e0b:	48 89 c7             	mov    %rax,%rdi
   11e0e:	e8 e5 fe ff ff       	call   11cf8 <_Z16update_particlesP10particle_tid>
   11e13:	90                   	nop
   11e14:	c9                   	leave  
   11e15:	c3                   	ret    

0000000000011e16 <_Z15alloc_particlesPP10particle_tiii>:
   11e16:	f3 0f 1e fa          	endbr64 
   11e1a:	55                   	push   %rbp
   11e1b:	48 89 e5             	mov    %rsp,%rbp
   11e1e:	48 83 ec 50          	sub    $0x50,%rsp
   11e22:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
   11e26:	89 75 c4             	mov    %esi,-0x3c(%rbp)
   11e29:	89 55 c0             	mov    %edx,-0x40(%rbp)
   11e2c:	89 4d bc             	mov    %ecx,-0x44(%rbp)
   11e2f:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   11e33:	48 8b 00             	mov    (%rax),%rax
   11e36:	48 85 c0             	test   %rax,%rax
   11e39:	74 28                	je     11e63 <_Z15alloc_particlesPP10particle_tiii+0x4d>
   11e3b:	48 8d 05 46 4b 00 00 	lea    0x4b46(%rip),%rax        # 16988 <_ZStL19piecewise_construct+0x208>
   11e42:	48 89 c1             	mov    %rax,%rcx
   11e45:	ba f6 00 00 00       	mov    $0xf6,%edx
   11e4a:	48 8d 05 6f 49 00 00 	lea    0x496f(%rip),%rax        # 167c0 <_ZStL19piecewise_construct+0x40>
   11e51:	48 89 c6             	mov    %rax,%rsi
   11e54:	48 8d 05 5f 4b 00 00 	lea    0x4b5f(%rip),%rax        # 169ba <_ZStL19piecewise_construct+0x23a>
   11e5b:	48 89 c7             	mov    %rax,%rdi
   11e5e:	e8 cd 82 ff ff       	call   a130 <__assert_fail@plt>
   11e63:	8b 45 c4             	mov    -0x3c(%rbp),%eax
   11e66:	48 98                	cltq   
   11e68:	48 c1 e0 06          	shl    $0x6,%rax
   11e6c:	48 89 c7             	mov    %rax,%rdi
   11e6f:	e8 ac 87 ff ff       	call   a620 <malloc@plt>
   11e74:	48 89 c2             	mov    %rax,%rdx
   11e77:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   11e7b:	48 89 10             	mov    %rdx,(%rax)
   11e7e:	8b 45 c4             	mov    -0x3c(%rbp),%eax
   11e81:	0f af 45 c0          	imul   -0x40(%rbp),%eax
   11e85:	48 98                	cltq   
   11e87:	48 c1 e0 03          	shl    $0x3,%rax
   11e8b:	48 89 c7             	mov    %rax,%rdi
   11e8e:	e8 8d 87 ff ff       	call   a620 <malloc@plt>
   11e93:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   11e97:	8b 45 c4             	mov    -0x3c(%rbp),%eax
   11e9a:	0f af 45 c0          	imul   -0x40(%rbp),%eax
   11e9e:	48 98                	cltq   
   11ea0:	48 c1 e0 03          	shl    $0x3,%rax
   11ea4:	48 89 c7             	mov    %rax,%rdi
   11ea7:	e8 74 87 ff ff       	call   a620 <malloc@plt>
   11eac:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   11eb0:	8b 45 c4             	mov    -0x3c(%rbp),%eax
   11eb3:	0f af 45 c0          	imul   -0x40(%rbp),%eax
   11eb7:	48 98                	cltq   
   11eb9:	48 c1 e0 03          	shl    $0x3,%rax
   11ebd:	48 89 c7             	mov    %rax,%rdi
   11ec0:	e8 5b 87 ff ff       	call   a620 <malloc@plt>
   11ec5:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   11ec9:	8b 45 c4             	mov    -0x3c(%rbp),%eax
   11ecc:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   11ed0:	48 98                	cltq   
   11ed2:	48 c1 e0 03          	shl    $0x3,%rax
   11ed6:	48 89 c7             	mov    %rax,%rdi
   11ed9:	e8 42 87 ff ff       	call   a620 <malloc@plt>
   11ede:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   11ee2:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   11ee6:	48 8b 00             	mov    (%rax),%rax
   11ee9:	48 85 c0             	test   %rax,%rax
   11eec:	75 28                	jne    11f16 <_Z15alloc_particlesPP10particle_tiii+0x100>
   11eee:	48 8d 05 93 4a 00 00 	lea    0x4a93(%rip),%rax        # 16988 <_ZStL19piecewise_construct+0x208>
   11ef5:	48 89 c1             	mov    %rax,%rcx
   11ef8:	ba fc 00 00 00       	mov    $0xfc,%edx
   11efd:	48 8d 05 bc 48 00 00 	lea    0x48bc(%rip),%rax        # 167c0 <_ZStL19piecewise_construct+0x40>
   11f04:	48 89 c6             	mov    %rax,%rsi
   11f07:	48 8d 05 c2 4a 00 00 	lea    0x4ac2(%rip),%rax        # 169d0 <_ZStL19piecewise_construct+0x250>
   11f0e:	48 89 c7             	mov    %rax,%rdi
   11f11:	e8 1a 82 ff ff       	call   a130 <__assert_fail@plt>
   11f16:	48 83 7d d8 00       	cmpq   $0x0,-0x28(%rbp)
   11f1b:	75 28                	jne    11f45 <_Z15alloc_particlesPP10particle_tiii+0x12f>
   11f1d:	48 8d 05 64 4a 00 00 	lea    0x4a64(%rip),%rax        # 16988 <_ZStL19piecewise_construct+0x208>
   11f24:	48 89 c1             	mov    %rax,%rcx
   11f27:	ba fd 00 00 00       	mov    $0xfd,%edx
   11f2c:	48 8d 05 8d 48 00 00 	lea    0x488d(%rip),%rax        # 167c0 <_ZStL19piecewise_construct+0x40>
   11f33:	48 89 c6             	mov    %rax,%rsi
   11f36:	48 8d 05 78 49 00 00 	lea    0x4978(%rip),%rax        # 168b5 <_ZStL19piecewise_construct+0x135>
   11f3d:	48 89 c7             	mov    %rax,%rdi
   11f40:	e8 eb 81 ff ff       	call   a130 <__assert_fail@plt>
   11f45:	48 83 7d e0 00       	cmpq   $0x0,-0x20(%rbp)
   11f4a:	75 28                	jne    11f74 <_Z15alloc_particlesPP10particle_tiii+0x15e>
   11f4c:	48 8d 05 35 4a 00 00 	lea    0x4a35(%rip),%rax        # 16988 <_ZStL19piecewise_construct+0x208>
   11f53:	48 89 c1             	mov    %rax,%rcx
   11f56:	ba fd 00 00 00       	mov    $0xfd,%edx
   11f5b:	48 8d 05 5e 48 00 00 	lea    0x485e(%rip),%rax        # 167c0 <_ZStL19piecewise_construct+0x40>
   11f62:	48 89 c6             	mov    %rax,%rsi
   11f65:	48 8d 05 58 49 00 00 	lea    0x4958(%rip),%rax        # 168c4 <_ZStL19piecewise_construct+0x144>
   11f6c:	48 89 c7             	mov    %rax,%rdi
   11f6f:	e8 bc 81 ff ff       	call   a130 <__assert_fail@plt>
   11f74:	48 83 7d e8 00       	cmpq   $0x0,-0x18(%rbp)
   11f79:	75 28                	jne    11fa3 <_Z15alloc_particlesPP10particle_tiii+0x18d>
   11f7b:	48 8d 05 06 4a 00 00 	lea    0x4a06(%rip),%rax        # 16988 <_ZStL19piecewise_construct+0x208>
   11f82:	48 89 c1             	mov    %rax,%rcx
   11f85:	ba fe 00 00 00       	mov    $0xfe,%edx
   11f8a:	48 8d 05 2f 48 00 00 	lea    0x482f(%rip),%rax        # 167c0 <_ZStL19piecewise_construct+0x40>
   11f91:	48 89 c6             	mov    %rax,%rsi
   11f94:	48 8d 05 38 49 00 00 	lea    0x4938(%rip),%rax        # 168d3 <_ZStL19piecewise_construct+0x153>
   11f9b:	48 89 c7             	mov    %rax,%rdi
   11f9e:	e8 8d 81 ff ff       	call   a130 <__assert_fail@plt>
   11fa3:	48 83 7d f0 00       	cmpq   $0x0,-0x10(%rbp)
   11fa8:	75 28                	jne    11fd2 <_Z15alloc_particlesPP10particle_tiii+0x1bc>
   11faa:	48 8d 05 d7 49 00 00 	lea    0x49d7(%rip),%rax        # 16988 <_ZStL19piecewise_construct+0x208>
   11fb1:	48 89 c1             	mov    %rax,%rcx
   11fb4:	ba fe 00 00 00       	mov    $0xfe,%edx
   11fb9:	48 8d 05 00 48 00 00 	lea    0x4800(%rip),%rax        # 167c0 <_ZStL19piecewise_construct+0x40>
   11fc0:	48 89 c6             	mov    %rax,%rsi
   11fc3:	48 8d 05 1c 4a 00 00 	lea    0x4a1c(%rip),%rax        # 169e6 <_ZStL19piecewise_construct+0x266>
   11fca:	48 89 c7             	mov    %rax,%rdi
   11fcd:	e8 5e 81 ff ff       	call   a130 <__assert_fail@plt>
   11fd2:	c7 45 d4 00 00 00 00 	movl   $0x0,-0x2c(%rbp)
   11fd9:	e9 ce 00 00 00       	jmp    120ac <_Z15alloc_particlesPP10particle_tiii+0x296>
   11fde:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   11fe2:	48 8b 10             	mov    (%rax),%rdx
   11fe5:	8b 45 d4             	mov    -0x2c(%rbp),%eax
   11fe8:	48 98                	cltq   
   11fea:	48 c1 e0 06          	shl    $0x6,%rax
   11fee:	48 01 d0             	add    %rdx,%rax
   11ff1:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   11ff5:	8b 45 d4             	mov    -0x2c(%rbp),%eax
   11ff8:	0f af 45 c0          	imul   -0x40(%rbp),%eax
   11ffc:	48 98                	cltq   
   11ffe:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12005:	00 
   12006:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   1200a:	48 01 c2             	add    %rax,%rdx
   1200d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12011:	48 89 50 08          	mov    %rdx,0x8(%rax)
   12015:	8b 45 d4             	mov    -0x2c(%rbp),%eax
   12018:	0f af 45 c0          	imul   -0x40(%rbp),%eax
   1201c:	48 98                	cltq   
   1201e:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12025:	00 
   12026:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   1202a:	48 01 c2             	add    %rax,%rdx
   1202d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12031:	48 89 50 10          	mov    %rdx,0x10(%rax)
   12035:	8b 45 d4             	mov    -0x2c(%rbp),%eax
   12038:	0f af 45 c0          	imul   -0x40(%rbp),%eax
   1203c:	48 98                	cltq   
   1203e:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12045:	00 
   12046:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1204a:	48 01 c2             	add    %rax,%rdx
   1204d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12051:	48 89 50 18          	mov    %rdx,0x18(%rax)
   12055:	8b 45 d4             	mov    -0x2c(%rbp),%eax
   12058:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   1205c:	48 98                	cltq   
   1205e:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12065:	00 
   12066:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   1206a:	48 01 c2             	add    %rax,%rdx
   1206d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12071:	48 89 50 20          	mov    %rdx,0x20(%rax)
   12075:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   12079:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1207d:	48 89 10             	mov    %rdx,(%rax)
   12080:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12084:	c6 40 38 00          	movb   $0x0,0x38(%rax)
   12088:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1208c:	8b 55 c0             	mov    -0x40(%rbp),%edx
   1208f:	89 50 30             	mov    %edx,0x30(%rax)
   12092:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12096:	8b 55 bc             	mov    -0x44(%rbp),%edx
   12099:	89 50 34             	mov    %edx,0x34(%rax)
   1209c:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   120a0:	48 c7 40 28 00 00 00 	movq   $0x0,0x28(%rax)
   120a7:	00 
   120a8:	83 45 d4 01          	addl   $0x1,-0x2c(%rbp)
   120ac:	8b 45 d4             	mov    -0x2c(%rbp),%eax
   120af:	3b 45 c4             	cmp    -0x3c(%rbp),%eax
   120b2:	0f 8c 26 ff ff ff    	jl     11fde <_Z15alloc_particlesPP10particle_tiii+0x1c8>
   120b8:	90                   	nop
   120b9:	90                   	nop
   120ba:	c9                   	leave  
   120bb:	c3                   	ret    

00000000000120bc <_Z15alloc_particlesP17chunk_particles_tiii>:
   120bc:	f3 0f 1e fa          	endbr64 
   120c0:	55                   	push   %rbp
   120c1:	48 89 e5             	mov    %rsp,%rbp
   120c4:	48 83 ec 30          	sub    $0x30,%rsp
   120c8:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   120cc:	89 75 e4             	mov    %esi,-0x1c(%rbp)
   120cf:	89 55 e0             	mov    %edx,-0x20(%rbp)
   120d2:	89 4d dc             	mov    %ecx,-0x24(%rbp)
   120d5:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   120dc:	00 00 
   120de:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   120e2:	31 c0                	xor    %eax,%eax
   120e4:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   120e8:	48 8b 00             	mov    (%rax),%rax
   120eb:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   120ef:	8b 4d dc             	mov    -0x24(%rbp),%ecx
   120f2:	8b 55 e0             	mov    -0x20(%rbp),%edx
   120f5:	8b 75 e4             	mov    -0x1c(%rbp),%esi
   120f8:	48 8d 45 f0          	lea    -0x10(%rbp),%rax
   120fc:	48 89 c7             	mov    %rax,%rdi
   120ff:	e8 12 fd ff ff       	call   11e16 <_Z15alloc_particlesPP10particle_tiii>
   12104:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   12108:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1210c:	48 89 10             	mov    %rdx,(%rax)
   1210f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12113:	8b 55 e4             	mov    -0x1c(%rbp),%edx
   12116:	89 50 08             	mov    %edx,0x8(%rax)
   12119:	90                   	nop
   1211a:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1211e:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   12125:	00 00 
   12127:	74 05                	je     1212e <_Z15alloc_particlesP17chunk_particles_tiii+0x72>
   12129:	e8 a2 83 ff ff       	call   a4d0 <__stack_chk_fail@plt>
   1212e:	c9                   	leave  
   1212f:	c3                   	ret    

0000000000012130 <_Z41__static_initialization_and_destruction_0ii>:
   12130:	f3 0f 1e fa          	endbr64 
   12134:	55                   	push   %rbp
   12135:	48 89 e5             	mov    %rsp,%rbp
   12138:	48 83 ec 10          	sub    $0x10,%rsp
   1213c:	89 7d fc             	mov    %edi,-0x4(%rbp)
   1213f:	89 75 f8             	mov    %esi,-0x8(%rbp)
   12142:	83 7d fc 01          	cmpl   $0x1,-0x4(%rbp)
   12146:	75 3b                	jne    12183 <_Z41__static_initialization_and_destruction_0ii+0x53>
   12148:	81 7d f8 ff ff 00 00 	cmpl   $0xffff,-0x8(%rbp)
   1214f:	75 32                	jne    12183 <_Z41__static_initialization_and_destruction_0ii+0x53>
   12151:	48 8d 05 d2 c0 00 00 	lea    0xc0d2(%rip),%rax        # 1e22a <_ZStL8__ioinit>
   12158:	48 89 c7             	mov    %rax,%rdi
   1215b:	e8 50 86 ff ff       	call   a7b0 <_ZNSt8ios_base4InitC1Ev@plt>
   12160:	48 8d 05 a1 ae 00 00 	lea    0xaea1(%rip),%rax        # 1d008 <__dso_handle>
   12167:	48 89 c2             	mov    %rax,%rdx
   1216a:	48 8d 05 b9 c0 00 00 	lea    0xc0b9(%rip),%rax        # 1e22a <_ZStL8__ioinit>
   12171:	48 89 c6             	mov    %rax,%rsi
   12174:	48 8b 05 7d ae 00 00 	mov    0xae7d(%rip),%rax        # 1cff8 <_ZNSt8ios_base4InitD1Ev@GLIBCXX_3.4>
   1217b:	48 89 c7             	mov    %rax,%rdi
   1217e:	e8 9d 81 ff ff       	call   a320 <__cxa_atexit@plt>
   12183:	90                   	nop
   12184:	c9                   	leave  
   12185:	c3                   	ret    

0000000000012186 <_GLOBAL__sub_I__Z14updateVelocityP10particle_td>:
   12186:	f3 0f 1e fa          	endbr64 
   1218a:	55                   	push   %rbp
   1218b:	48 89 e5             	mov    %rsp,%rbp
   1218e:	be ff ff 00 00       	mov    $0xffff,%esi
   12193:	bf 01 00 00 00       	mov    $0x1,%edi
   12198:	e8 93 ff ff ff       	call   12130 <_Z41__static_initialization_and_destruction_0ii>
   1219d:	5d                   	pop    %rbp
   1219e:	c3                   	ret    

000000000001219f <_ZStanSt13_Ios_FmtflagsS_>:
   1219f:	f3 0f 1e fa          	endbr64 
   121a3:	55                   	push   %rbp
   121a4:	48 89 e5             	mov    %rsp,%rbp
   121a7:	89 7d fc             	mov    %edi,-0x4(%rbp)
   121aa:	89 75 f8             	mov    %esi,-0x8(%rbp)
   121ad:	8b 45 fc             	mov    -0x4(%rbp),%eax
   121b0:	23 45 f8             	and    -0x8(%rbp),%eax
   121b3:	5d                   	pop    %rbp
   121b4:	c3                   	ret    

00000000000121b5 <_ZStorSt13_Ios_FmtflagsS_>:
   121b5:	f3 0f 1e fa          	endbr64 
   121b9:	55                   	push   %rbp
   121ba:	48 89 e5             	mov    %rsp,%rbp
   121bd:	89 7d fc             	mov    %edi,-0x4(%rbp)
   121c0:	89 75 f8             	mov    %esi,-0x8(%rbp)
   121c3:	8b 45 fc             	mov    -0x4(%rbp),%eax
   121c6:	0b 45 f8             	or     -0x8(%rbp),%eax
   121c9:	5d                   	pop    %rbp
   121ca:	c3                   	ret    

00000000000121cb <_ZStcoSt13_Ios_Fmtflags>:
   121cb:	f3 0f 1e fa          	endbr64 
   121cf:	55                   	push   %rbp
   121d0:	48 89 e5             	mov    %rsp,%rbp
   121d3:	89 7d fc             	mov    %edi,-0x4(%rbp)
   121d6:	8b 45 fc             	mov    -0x4(%rbp),%eax
   121d9:	f7 d0                	not    %eax
   121db:	5d                   	pop    %rbp
   121dc:	c3                   	ret    

00000000000121dd <_ZStoRRSt13_Ios_FmtflagsS_>:
   121dd:	f3 0f 1e fa          	endbr64 
   121e1:	55                   	push   %rbp
   121e2:	48 89 e5             	mov    %rsp,%rbp
   121e5:	48 83 ec 10          	sub    $0x10,%rsp
   121e9:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   121ed:	89 75 f4             	mov    %esi,-0xc(%rbp)
   121f0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   121f4:	8b 00                	mov    (%rax),%eax
   121f6:	8b 55 f4             	mov    -0xc(%rbp),%edx
   121f9:	89 d6                	mov    %edx,%esi
   121fb:	89 c7                	mov    %eax,%edi
   121fd:	e8 b3 ff ff ff       	call   121b5 <_ZStorSt13_Ios_FmtflagsS_>
   12202:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   12206:	89 02                	mov    %eax,(%rdx)
   12208:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1220c:	c9                   	leave  
   1220d:	c3                   	ret    

000000000001220e <_ZStaNRSt13_Ios_FmtflagsS_>:
   1220e:	f3 0f 1e fa          	endbr64 
   12212:	55                   	push   %rbp
   12213:	48 89 e5             	mov    %rsp,%rbp
   12216:	48 83 ec 10          	sub    $0x10,%rsp
   1221a:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   1221e:	89 75 f4             	mov    %esi,-0xc(%rbp)
   12221:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12225:	8b 00                	mov    (%rax),%eax
   12227:	8b 55 f4             	mov    -0xc(%rbp),%edx
   1222a:	89 d6                	mov    %edx,%esi
   1222c:	89 c7                	mov    %eax,%edi
   1222e:	e8 6c ff ff ff       	call   1219f <_ZStanSt13_Ios_FmtflagsS_>
   12233:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   12237:	89 02                	mov    %eax,(%rdx)
   12239:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1223d:	c9                   	leave  
   1223e:	c3                   	ret    

000000000001223f <_ZStorSt13_Ios_OpenmodeS_>:
   1223f:	f3 0f 1e fa          	endbr64 
   12243:	55                   	push   %rbp
   12244:	48 89 e5             	mov    %rsp,%rbp
   12247:	89 7d fc             	mov    %edi,-0x4(%rbp)
   1224a:	89 75 f8             	mov    %esi,-0x8(%rbp)
   1224d:	8b 45 fc             	mov    -0x4(%rbp),%eax
   12250:	0b 45 f8             	or     -0x8(%rbp),%eax
   12253:	5d                   	pop    %rbp
   12254:	c3                   	ret    
   12255:	90                   	nop

0000000000012256 <_ZNSt8ios_base4setfESt13_Ios_FmtflagsS0_>:
   12256:	f3 0f 1e fa          	endbr64 
   1225a:	55                   	push   %rbp
   1225b:	48 89 e5             	mov    %rsp,%rbp
   1225e:	48 83 ec 20          	sub    $0x20,%rsp
   12262:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   12266:	89 75 e4             	mov    %esi,-0x1c(%rbp)
   12269:	89 55 e0             	mov    %edx,-0x20(%rbp)
   1226c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12270:	8b 40 18             	mov    0x18(%rax),%eax
   12273:	89 45 fc             	mov    %eax,-0x4(%rbp)
   12276:	8b 45 e0             	mov    -0x20(%rbp),%eax
   12279:	89 c7                	mov    %eax,%edi
   1227b:	e8 4b ff ff ff       	call   121cb <_ZStcoSt13_Ios_Fmtflags>
   12280:	89 c2                	mov    %eax,%edx
   12282:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12286:	48 83 c0 18          	add    $0x18,%rax
   1228a:	89 d6                	mov    %edx,%esi
   1228c:	48 89 c7             	mov    %rax,%rdi
   1228f:	e8 7a ff ff ff       	call   1220e <_ZStaNRSt13_Ios_FmtflagsS_>
   12294:	8b 55 e0             	mov    -0x20(%rbp),%edx
   12297:	8b 45 e4             	mov    -0x1c(%rbp),%eax
   1229a:	89 d6                	mov    %edx,%esi
   1229c:	89 c7                	mov    %eax,%edi
   1229e:	e8 fc fe ff ff       	call   1219f <_ZStanSt13_Ios_FmtflagsS_>
   122a3:	89 c2                	mov    %eax,%edx
   122a5:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   122a9:	48 83 c0 18          	add    $0x18,%rax
   122ad:	89 d6                	mov    %edx,%esi
   122af:	48 89 c7             	mov    %rax,%rdi
   122b2:	e8 26 ff ff ff       	call   121dd <_ZStoRRSt13_Ios_FmtflagsS_>
   122b7:	8b 45 fc             	mov    -0x4(%rbp),%eax
   122ba:	c9                   	leave  
   122bb:	c3                   	ret    

00000000000122bc <_ZSt10scientificRSt8ios_base>:
   122bc:	f3 0f 1e fa          	endbr64 
   122c0:	55                   	push   %rbp
   122c1:	48 89 e5             	mov    %rsp,%rbp
   122c4:	48 83 ec 10          	sub    $0x10,%rsp
   122c8:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   122cc:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   122d0:	ba 04 01 00 00       	mov    $0x104,%edx
   122d5:	be 00 01 00 00       	mov    $0x100,%esi
   122da:	48 89 c7             	mov    %rax,%rdi
   122dd:	e8 74 ff ff ff       	call   12256 <_ZNSt8ios_base4setfESt13_Ios_FmtflagsS0_>
   122e2:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   122e6:	c9                   	leave  
   122e7:	c3                   	ret    

00000000000122e8 <_Z25gravityUpdateAccelerationP10particle_tS0_>:
   122e8:	f3 0f 1e fa          	endbr64 
   122ec:	55                   	push   %rbp
   122ed:	48 89 e5             	mov    %rsp,%rbp
   122f0:	48 83 ec 70          	sub    $0x70,%rsp
   122f4:	48 89 7d 98          	mov    %rdi,-0x68(%rbp)
   122f8:	48 89 75 90          	mov    %rsi,-0x70(%rbp)
   122fc:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   12303:	00 00 
   12305:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   12309:	31 c0                	xor    %eax,%eax
   1230b:	f2 0f 10 05 f5 46 00 	movsd  0x46f5(%rip),%xmm0        # 16a08 <_ZStL19piecewise_construct+0x8>
   12312:	00 
   12313:	f2 0f 11 45 c8       	movsd  %xmm0,-0x38(%rbp)
   12318:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1231c:	48 8b 40 20          	mov    0x20(%rax),%rax
   12320:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   12324:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   12328:	48 8b 40 20          	mov    0x20(%rax),%rax
   1232c:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   12330:	48 c7 45 b0 00 00 00 	movq   $0x0,-0x50(%rbp)
   12337:	00 
   12338:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1233c:	8b 50 30             	mov    0x30(%rax),%edx
   1233f:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   12343:	48 8b 70 08          	mov    0x8(%rax),%rsi
   12347:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1234b:	48 8b 40 08          	mov    0x8(%rax),%rax
   1234f:	48 8d 4d b0          	lea    -0x50(%rbp),%rcx
   12353:	48 89 c7             	mov    %rax,%rdi
   12356:	e8 21 2a 00 00       	call   14d7c <_Z16get_distance_vecPKdS0_iPPd>
   1235b:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1235f:	8b 50 30             	mov    0x30(%rax),%edx
   12362:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   12366:	89 d6                	mov    %edx,%esi
   12368:	48 89 c7             	mov    %rax,%rdi
   1236b:	e8 a4 2a 00 00       	call   14e14 <_Z10get_lengthPKdi>
   12370:	66 48 0f 7e c0       	movq   %xmm0,%rax
   12375:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   12379:	f2 0f 10 45 e0       	movsd  -0x20(%rbp),%xmm0
   1237e:	f2 0f 59 c0          	mulsd  %xmm0,%xmm0
   12382:	f2 0f 10 4d e0       	movsd  -0x20(%rbp),%xmm1
   12387:	f2 0f 59 c1          	mulsd  %xmm1,%xmm0
   1238b:	f2 0f 11 45 e0       	movsd  %xmm0,-0x20(%rbp)
   12390:	f2 0f 10 05 70 46 00 	movsd  0x4670(%rip),%xmm0        # 16a08 <_ZStL19piecewise_construct+0x8>
   12397:	00 
   12398:	66 0f 28 c8          	movapd %xmm0,%xmm1
   1239c:	f2 0f 5e 4d e0       	divsd  -0x20(%rbp),%xmm1
   123a1:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   123a5:	f2 0f 10 00          	movsd  (%rax),%xmm0
   123a9:	f2 0f 59 c1          	mulsd  %xmm1,%xmm0
   123ad:	f2 0f 11 45 e8       	movsd  %xmm0,-0x18(%rbp)
   123b2:	f2 0f 10 05 4e 46 00 	movsd  0x464e(%rip),%xmm0        # 16a08 <_ZStL19piecewise_construct+0x8>
   123b9:	00 
   123ba:	66 0f 28 c8          	movapd %xmm0,%xmm1
   123be:	f2 0f 5e 4d e0       	divsd  -0x20(%rbp),%xmm1
   123c3:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   123c7:	f2 0f 10 00          	movsd  (%rax),%xmm0
   123cb:	f2 0f 59 c1          	mulsd  %xmm1,%xmm0
   123cf:	f2 0f 11 45 f0       	movsd  %xmm0,-0x10(%rbp)
   123d4:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   123d8:	48 8b 40 18          	mov    0x18(%rax),%rax
   123dc:	48 89 45 b8          	mov    %rax,-0x48(%rbp)
   123e0:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   123e4:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   123e8:	c7 45 a8 00 00 00 00 	movl   $0x0,-0x58(%rbp)
   123ef:	eb 35                	jmp    12426 <_Z25gravityUpdateAccelerationP10particle_tS0_+0x13e>
   123f1:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   123f5:	48 8d 50 08          	lea    0x8(%rax),%rdx
   123f9:	48 89 55 c0          	mov    %rdx,-0x40(%rbp)
   123fd:	f2 0f 10 00          	movsd  (%rax),%xmm0
   12401:	66 0f 28 c8          	movapd %xmm0,%xmm1
   12405:	f2 0f 59 4d e8       	mulsd  -0x18(%rbp),%xmm1
   1240a:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   1240e:	48 8d 50 08          	lea    0x8(%rax),%rdx
   12412:	48 89 55 b8          	mov    %rdx,-0x48(%rbp)
   12416:	f2 0f 10 00          	movsd  (%rax),%xmm0
   1241a:	f2 0f 58 c1          	addsd  %xmm1,%xmm0
   1241e:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   12422:	83 45 a8 01          	addl   $0x1,-0x58(%rbp)
   12426:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1242a:	8b 40 30             	mov    0x30(%rax),%eax
   1242d:	39 45 a8             	cmp    %eax,-0x58(%rbp)
   12430:	7c bf                	jl     123f1 <_Z25gravityUpdateAccelerationP10particle_tS0_+0x109>
   12432:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   12436:	48 8b 40 18          	mov    0x18(%rax),%rax
   1243a:	48 89 45 b8          	mov    %rax,-0x48(%rbp)
   1243e:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   12442:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   12446:	c7 45 ac 00 00 00 00 	movl   $0x0,-0x54(%rbp)
   1244d:	eb 35                	jmp    12484 <_Z25gravityUpdateAccelerationP10particle_tS0_+0x19c>
   1244f:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   12453:	48 8d 50 08          	lea    0x8(%rax),%rdx
   12457:	48 89 55 c0          	mov    %rdx,-0x40(%rbp)
   1245b:	f2 0f 10 00          	movsd  (%rax),%xmm0
   1245f:	66 0f 28 c8          	movapd %xmm0,%xmm1
   12463:	f2 0f 59 4d f0       	mulsd  -0x10(%rbp),%xmm1
   12468:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   1246c:	48 8d 50 08          	lea    0x8(%rax),%rdx
   12470:	48 89 55 b8          	mov    %rdx,-0x48(%rbp)
   12474:	f2 0f 10 00          	movsd  (%rax),%xmm0
   12478:	f2 0f 5c c1          	subsd  %xmm1,%xmm0
   1247c:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   12480:	83 45 ac 01          	addl   $0x1,-0x54(%rbp)
   12484:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   12488:	8b 40 30             	mov    0x30(%rax),%eax
   1248b:	39 45 ac             	cmp    %eax,-0x54(%rbp)
   1248e:	7c bf                	jl     1244f <_Z25gravityUpdateAccelerationP10particle_tS0_+0x167>
   12490:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   12494:	48 89 c7             	mov    %rax,%rdi
   12497:	e8 54 80 ff ff       	call   a4f0 <free@plt>
   1249c:	90                   	nop
   1249d:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   124a1:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   124a8:	00 00 
   124aa:	74 05                	je     124b1 <_Z25gravityUpdateAccelerationP10particle_tS0_+0x1c9>
   124ac:	e8 1f 80 ff ff       	call   a4d0 <__stack_chk_fail@plt>
   124b1:	c9                   	leave  
   124b2:	c3                   	ret    

00000000000124b3 <_Z18setupCommunicatorsP10topology_tii>:
   124b3:	f3 0f 1e fa          	endbr64 
   124b7:	55                   	push   %rbp
   124b8:	48 89 e5             	mov    %rsp,%rbp
   124bb:	48 83 ec 20          	sub    $0x20,%rsp
   124bf:	48 89 7d e8          	mov    %rdi,-0x18(%rbp)
   124c3:	89 75 e4             	mov    %esi,-0x1c(%rbp)
   124c6:	89 55 e0             	mov    %edx,-0x20(%rbp)
   124c9:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   124d0:	00 00 
   124d2:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   124d6:	31 c0                	xor    %eax,%eax
   124d8:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   124dc:	8b 55 e4             	mov    -0x1c(%rbp),%edx
   124df:	89 50 08             	mov    %edx,0x8(%rax)
   124e2:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   124e6:	8b 55 e0             	mov    -0x20(%rbp),%edx
   124e9:	89 50 0c             	mov    %edx,0xc(%rax)
   124ec:	48 8d 45 f0          	lea    -0x10(%rbp),%rax
   124f0:	48 89 c6             	mov    %rax,%rsi
   124f3:	48 8d 05 26 b5 00 00 	lea    0xb526(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   124fa:	48 89 c7             	mov    %rax,%rdi
   124fd:	e8 de 7f ff ff       	call   a4e0 <MPI_Comm_size@plt>
   12502:	8b 45 e4             	mov    -0x1c(%rbp),%eax
   12505:	0f af 45 e0          	imul   -0x20(%rbp),%eax
   12509:	89 c2                	mov    %eax,%edx
   1250b:	8b 45 f0             	mov    -0x10(%rbp),%eax
   1250e:	39 c2                	cmp    %eax,%edx
   12510:	74 28                	je     1253a <_Z18setupCommunicatorsP10topology_tii+0x87>
   12512:	48 8d 05 37 47 00 00 	lea    0x4737(%rip),%rax        # 16c50 <_ZN3MPIL11LOCK_SHAREDE+0x8>
   12519:	48 89 c1             	mov    %rax,%rcx
   1251c:	ba 0f 00 00 00       	mov    $0xf,%edx
   12521:	48 8d 05 58 47 00 00 	lea    0x4758(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12528:	48 89 c6             	mov    %rax,%rsi
   1252b:	48 8d 05 a6 47 00 00 	lea    0x47a6(%rip),%rax        # 16cd8 <_ZN3MPIL11LOCK_SHAREDE+0x90>
   12532:	48 89 c7             	mov    %rax,%rdi
   12535:	e8 f6 7b ff ff       	call   a130 <__assert_fail@plt>
   1253a:	48 8d 45 f4          	lea    -0xc(%rbp),%rax
   1253e:	48 89 c6             	mov    %rax,%rsi
   12541:	48 8d 05 d8 b4 00 00 	lea    0xb4d8(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   12548:	48 89 c7             	mov    %rax,%rdi
   1254b:	e8 a0 81 ff ff       	call   a6f0 <MPI_Comm_rank@plt>
   12550:	8b 45 f4             	mov    -0xc(%rbp),%eax
   12553:	99                   	cltd   
   12554:	f7 7d e4             	idivl  -0x1c(%rbp)
   12557:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1255b:	89 10                	mov    %edx,(%rax)
   1255d:	8b 45 f4             	mov    -0xc(%rbp),%eax
   12560:	99                   	cltd   
   12561:	f7 7d e4             	idivl  -0x1c(%rbp)
   12564:	89 c2                	mov    %eax,%edx
   12566:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1256a:	89 50 04             	mov    %edx,0x4(%rax)
   1256d:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12571:	8b 40 04             	mov    0x4(%rax),%eax
   12574:	85 c0                	test   %eax,%eax
   12576:	75 28                	jne    125a0 <_Z18setupCommunicatorsP10topology_tii+0xed>
   12578:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1257c:	48 8d 48 20          	lea    0x20(%rax),%rcx
   12580:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12584:	8b 10                	mov    (%rax),%edx
   12586:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1258a:	8b 40 04             	mov    0x4(%rax),%eax
   1258d:	89 c6                	mov    %eax,%esi
   1258f:	48 8d 05 8a b4 00 00 	lea    0xb48a(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   12596:	48 89 c7             	mov    %rax,%rdi
   12599:	e8 c2 82 ff ff       	call   a860 <MPI_Comm_split@plt>
   1259e:	eb 27                	jmp    125c7 <_Z18setupCommunicatorsP10topology_tii+0x114>
   125a0:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   125a4:	48 8d 50 20          	lea    0x20(%rax),%rdx
   125a8:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   125ac:	8b 00                	mov    (%rax),%eax
   125ae:	48 89 d1             	mov    %rdx,%rcx
   125b1:	89 c2                	mov    %eax,%edx
   125b3:	be 02 80 ff ff       	mov    $0xffff8002,%esi
   125b8:	48 8d 05 61 b4 00 00 	lea    0xb461(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   125bf:	48 89 c7             	mov    %rax,%rdi
   125c2:	e8 99 82 ff ff       	call   a860 <MPI_Comm_split@plt>
   125c7:	48 8d 05 52 b4 00 00 	lea    0xb452(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   125ce:	48 89 c7             	mov    %rax,%rdi
   125d1:	e8 ea 7b ff ff       	call   a1c0 <MPI_Barrier@plt>
   125d6:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   125da:	8b 00                	mov    (%rax),%eax
   125dc:	85 c0                	test   %eax,%eax
   125de:	75 28                	jne    12608 <_Z18setupCommunicatorsP10topology_tii+0x155>
   125e0:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   125e4:	48 8d 48 28          	lea    0x28(%rax),%rcx
   125e8:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   125ec:	8b 50 04             	mov    0x4(%rax),%edx
   125ef:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   125f3:	8b 00                	mov    (%rax),%eax
   125f5:	89 c6                	mov    %eax,%esi
   125f7:	48 8d 05 22 b4 00 00 	lea    0xb422(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   125fe:	48 89 c7             	mov    %rax,%rdi
   12601:	e8 5a 82 ff ff       	call   a860 <MPI_Comm_split@plt>
   12606:	eb 28                	jmp    12630 <_Z18setupCommunicatorsP10topology_tii+0x17d>
   12608:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1260c:	48 8d 50 28          	lea    0x28(%rax),%rdx
   12610:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12614:	8b 40 04             	mov    0x4(%rax),%eax
   12617:	48 89 d1             	mov    %rdx,%rcx
   1261a:	89 c2                	mov    %eax,%edx
   1261c:	be 02 80 ff ff       	mov    $0xffff8002,%esi
   12621:	48 8d 05 f8 b3 00 00 	lea    0xb3f8(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   12628:	48 89 c7             	mov    %rax,%rdi
   1262b:	e8 30 82 ff ff       	call   a860 <MPI_Comm_split@plt>
   12630:	48 8d 05 e9 b3 00 00 	lea    0xb3e9(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   12637:	48 89 c7             	mov    %rax,%rdi
   1263a:	e8 81 7b ff ff       	call   a1c0 <MPI_Barrier@plt>
   1263f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12643:	48 8d 48 10          	lea    0x10(%rax),%rcx
   12647:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1264b:	8b 50 04             	mov    0x4(%rax),%edx
   1264e:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12652:	8b 00                	mov    (%rax),%eax
   12654:	89 c6                	mov    %eax,%esi
   12656:	48 8d 05 c3 b3 00 00 	lea    0xb3c3(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   1265d:	48 89 c7             	mov    %rax,%rdi
   12660:	e8 fb 81 ff ff       	call   a860 <MPI_Comm_split@plt>
   12665:	48 8d 05 b4 b3 00 00 	lea    0xb3b4(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   1266c:	48 89 c7             	mov    %rax,%rdi
   1266f:	e8 4c 7b ff ff       	call   a1c0 <MPI_Barrier@plt>
   12674:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12678:	48 8d 48 18          	lea    0x18(%rax),%rcx
   1267c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12680:	8b 10                	mov    (%rax),%edx
   12682:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12686:	8b 40 04             	mov    0x4(%rax),%eax
   12689:	89 c6                	mov    %eax,%esi
   1268b:	48 8d 05 8e b3 00 00 	lea    0xb38e(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   12692:	48 89 c7             	mov    %rax,%rdi
   12695:	e8 c6 81 ff ff       	call   a860 <MPI_Comm_split@plt>
   1269a:	48 8d 05 7f b3 00 00 	lea    0xb37f(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   126a1:	48 89 c7             	mov    %rax,%rdi
   126a4:	e8 17 7b ff ff       	call   a1c0 <MPI_Barrier@plt>
   126a9:	90                   	nop
   126aa:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   126ae:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   126b5:	00 00 
   126b7:	74 05                	je     126be <_Z18setupCommunicatorsP10topology_tii+0x20b>
   126b9:	e8 12 7e ff ff       	call   a4d0 <__stack_chk_fail@plt>
   126be:	c9                   	leave  
   126bf:	c3                   	ret    

00000000000126c0 <_Z7paddingPP10particle_tPii>:
   126c0:	f3 0f 1e fa          	endbr64 
   126c4:	55                   	push   %rbp
   126c5:	48 89 e5             	mov    %rsp,%rbp
   126c8:	48 81 ec a0 00 00 00 	sub    $0xa0,%rsp
   126cf:	48 89 bd 78 ff ff ff 	mov    %rdi,-0x88(%rbp)
   126d6:	48 89 b5 70 ff ff ff 	mov    %rsi,-0x90(%rbp)
   126dd:	89 95 6c ff ff ff    	mov    %edx,-0x94(%rbp)
   126e3:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   126ea:	8b 00                	mov    (%rax),%eax
   126ec:	89 45 a0             	mov    %eax,-0x60(%rbp)
   126ef:	8b 95 6c ff ff ff    	mov    -0x94(%rbp),%edx
   126f5:	8b 45 a0             	mov    -0x60(%rbp),%eax
   126f8:	89 d6                	mov    %edx,%esi
   126fa:	89 c7                	mov    %eax,%edi
   126fc:	e8 58 26 00 00       	call   14d59 <_Z16sizeAfterPaddingii>
   12701:	89 45 a4             	mov    %eax,-0x5c(%rbp)
   12704:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   1270b:	48 8b 00             	mov    (%rax),%rax
   1270e:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
   12712:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   12716:	48 8b 40 08          	mov    0x8(%rax),%rax
   1271a:	48 89 45 b8          	mov    %rax,-0x48(%rbp)
   1271e:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   12722:	48 8b 40 10          	mov    0x10(%rax),%rax
   12726:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   1272a:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   1272e:	48 8b 40 18          	mov    0x18(%rax),%rax
   12732:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   12736:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   1273a:	48 8b 40 20          	mov    0x20(%rax),%rax
   1273e:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   12742:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   12746:	8b 40 30             	mov    0x30(%rax),%eax
   12749:	89 45 a8             	mov    %eax,-0x58(%rbp)
   1274c:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   12750:	8b 40 34             	mov    0x34(%rax),%eax
   12753:	89 45 ac             	mov    %eax,-0x54(%rbp)
   12756:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   12759:	48 98                	cltq   
   1275b:	48 c1 e0 06          	shl    $0x6,%rax
   1275f:	48 89 c2             	mov    %rax,%rdx
   12762:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   12766:	48 89 d6             	mov    %rdx,%rsi
   12769:	48 89 c7             	mov    %rax,%rdi
   1276c:	e8 2f 81 ff ff       	call   a8a0 <realloc@plt>
   12771:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   12775:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   12778:	0f af 45 a8          	imul   -0x58(%rbp),%eax
   1277c:	48 98                	cltq   
   1277e:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12785:	00 
   12786:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   1278a:	48 89 d6             	mov    %rdx,%rsi
   1278d:	48 89 c7             	mov    %rax,%rdi
   12790:	e8 0b 81 ff ff       	call   a8a0 <realloc@plt>
   12795:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   12799:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   1279c:	0f af 45 a8          	imul   -0x58(%rbp),%eax
   127a0:	48 98                	cltq   
   127a2:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   127a9:	00 
   127aa:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   127ae:	48 89 d6             	mov    %rdx,%rsi
   127b1:	48 89 c7             	mov    %rax,%rdi
   127b4:	e8 e7 80 ff ff       	call   a8a0 <realloc@plt>
   127b9:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   127bd:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   127c0:	0f af 45 a8          	imul   -0x58(%rbp),%eax
   127c4:	48 98                	cltq   
   127c6:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   127cd:	00 
   127ce:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   127d2:	48 89 d6             	mov    %rdx,%rsi
   127d5:	48 89 c7             	mov    %rax,%rdi
   127d8:	e8 c3 80 ff ff       	call   a8a0 <realloc@plt>
   127dd:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   127e1:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   127e4:	0f af 45 ac          	imul   -0x54(%rbp),%eax
   127e8:	48 98                	cltq   
   127ea:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   127f1:	00 
   127f2:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   127f6:	48 89 d6             	mov    %rdx,%rsi
   127f9:	48 89 c7             	mov    %rax,%rdi
   127fc:	e8 9f 80 ff ff       	call   a8a0 <realloc@plt>
   12801:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   12805:	48 83 7d d8 00       	cmpq   $0x0,-0x28(%rbp)
   1280a:	75 28                	jne    12834 <_Z7paddingPP10particle_tPii+0x174>
   1280c:	48 8d 05 e5 44 00 00 	lea    0x44e5(%rip),%rax        # 16cf8 <_ZN3MPIL11LOCK_SHAREDE+0xb0>
   12813:	48 89 c1             	mov    %rax,%rcx
   12816:	ba 3c 00 00 00       	mov    $0x3c,%edx
   1281b:	48 8d 05 5e 44 00 00 	lea    0x445e(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12822:	48 89 c6             	mov    %rax,%rsi
   12825:	48 8d 05 f2 44 00 00 	lea    0x44f2(%rip),%rax        # 16d1e <_ZN3MPIL11LOCK_SHAREDE+0xd6>
   1282c:	48 89 c7             	mov    %rax,%rdi
   1282f:	e8 fc 78 ff ff       	call   a130 <__assert_fail@plt>
   12834:	48 83 7d e0 00       	cmpq   $0x0,-0x20(%rbp)
   12839:	75 28                	jne    12863 <_Z7paddingPP10particle_tPii+0x1a3>
   1283b:	48 8d 05 b6 44 00 00 	lea    0x44b6(%rip),%rax        # 16cf8 <_ZN3MPIL11LOCK_SHAREDE+0xb0>
   12842:	48 89 c1             	mov    %rax,%rcx
   12845:	ba 3c 00 00 00       	mov    $0x3c,%edx
   1284a:	48 8d 05 2f 44 00 00 	lea    0x442f(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12851:	48 89 c6             	mov    %rax,%rsi
   12854:	48 8d 05 d3 44 00 00 	lea    0x44d3(%rip),%rax        # 16d2e <_ZN3MPIL11LOCK_SHAREDE+0xe6>
   1285b:	48 89 c7             	mov    %rax,%rdi
   1285e:	e8 cd 78 ff ff       	call   a130 <__assert_fail@plt>
   12863:	48 83 7d e8 00       	cmpq   $0x0,-0x18(%rbp)
   12868:	75 28                	jne    12892 <_Z7paddingPP10particle_tPii+0x1d2>
   1286a:	48 8d 05 87 44 00 00 	lea    0x4487(%rip),%rax        # 16cf8 <_ZN3MPIL11LOCK_SHAREDE+0xb0>
   12871:	48 89 c1             	mov    %rax,%rcx
   12874:	ba 3d 00 00 00       	mov    $0x3d,%edx
   12879:	48 8d 05 00 44 00 00 	lea    0x4400(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12880:	48 89 c6             	mov    %rax,%rsi
   12883:	48 8d 05 b6 44 00 00 	lea    0x44b6(%rip),%rax        # 16d40 <_ZN3MPIL11LOCK_SHAREDE+0xf8>
   1288a:	48 89 c7             	mov    %rax,%rdi
   1288d:	e8 9e 78 ff ff       	call   a130 <__assert_fail@plt>
   12892:	48 83 7d f0 00       	cmpq   $0x0,-0x10(%rbp)
   12897:	75 28                	jne    128c1 <_Z7paddingPP10particle_tPii+0x201>
   12899:	48 8d 05 58 44 00 00 	lea    0x4458(%rip),%rax        # 16cf8 <_ZN3MPIL11LOCK_SHAREDE+0xb0>
   128a0:	48 89 c1             	mov    %rax,%rcx
   128a3:	ba 3d 00 00 00       	mov    $0x3d,%edx
   128a8:	48 8d 05 d1 43 00 00 	lea    0x43d1(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   128af:	48 89 c6             	mov    %rax,%rsi
   128b2:	48 8d 05 99 44 00 00 	lea    0x4499(%rip),%rax        # 16d52 <_ZN3MPIL11LOCK_SHAREDE+0x10a>
   128b9:	48 89 c7             	mov    %rax,%rdi
   128bc:	e8 6f 78 ff ff       	call   a130 <__assert_fail@plt>
   128c1:	8b 45 a0             	mov    -0x60(%rbp),%eax
   128c4:	89 45 8c             	mov    %eax,-0x74(%rbp)
   128c7:	e9 98 00 00 00       	jmp    12964 <_Z7paddingPP10particle_tPii+0x2a4>
   128cc:	c7 45 90 00 00 00 00 	movl   $0x0,-0x70(%rbp)
   128d3:	eb 7f                	jmp    12954 <_Z7paddingPP10particle_tPii+0x294>
   128d5:	8b 45 8c             	mov    -0x74(%rbp),%eax
   128d8:	0f af 45 a8          	imul   -0x58(%rbp),%eax
   128dc:	48 63 d0             	movslq %eax,%rdx
   128df:	8b 45 90             	mov    -0x70(%rbp),%eax
   128e2:	48 98                	cltq   
   128e4:	48 01 d0             	add    %rdx,%rax
   128e7:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   128ee:	00 
   128ef:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   128f3:	48 01 d0             	add    %rdx,%rax
   128f6:	66 0f ef c0          	pxor   %xmm0,%xmm0
   128fa:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   128fe:	8b 45 8c             	mov    -0x74(%rbp),%eax
   12901:	0f af 45 a8          	imul   -0x58(%rbp),%eax
   12905:	48 63 d0             	movslq %eax,%rdx
   12908:	8b 45 90             	mov    -0x70(%rbp),%eax
   1290b:	48 98                	cltq   
   1290d:	48 01 d0             	add    %rdx,%rax
   12910:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12917:	00 
   12918:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1291c:	48 01 d0             	add    %rdx,%rax
   1291f:	66 0f ef c0          	pxor   %xmm0,%xmm0
   12923:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   12927:	8b 45 8c             	mov    -0x74(%rbp),%eax
   1292a:	0f af 45 a8          	imul   -0x58(%rbp),%eax
   1292e:	48 63 d0             	movslq %eax,%rdx
   12931:	8b 45 90             	mov    -0x70(%rbp),%eax
   12934:	48 98                	cltq   
   12936:	48 01 d0             	add    %rdx,%rax
   12939:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12940:	00 
   12941:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   12945:	48 01 d0             	add    %rdx,%rax
   12948:	66 0f ef c0          	pxor   %xmm0,%xmm0
   1294c:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   12950:	83 45 90 01          	addl   $0x1,-0x70(%rbp)
   12954:	8b 45 90             	mov    -0x70(%rbp),%eax
   12957:	3b 45 a8             	cmp    -0x58(%rbp),%eax
   1295a:	0f 8c 75 ff ff ff    	jl     128d5 <_Z7paddingPP10particle_tPii+0x215>
   12960:	83 45 8c 01          	addl   $0x1,-0x74(%rbp)
   12964:	8b 45 8c             	mov    -0x74(%rbp),%eax
   12967:	3b 45 a4             	cmp    -0x5c(%rbp),%eax
   1296a:	0f 8c 5c ff ff ff    	jl     128cc <_Z7paddingPP10particle_tPii+0x20c>
   12970:	8b 45 a0             	mov    -0x60(%rbp),%eax
   12973:	89 45 94             	mov    %eax,-0x6c(%rbp)
   12976:	eb 42                	jmp    129ba <_Z7paddingPP10particle_tPii+0x2fa>
   12978:	c7 45 98 00 00 00 00 	movl   $0x0,-0x68(%rbp)
   1297f:	eb 2d                	jmp    129ae <_Z7paddingPP10particle_tPii+0x2ee>
   12981:	8b 45 94             	mov    -0x6c(%rbp),%eax
   12984:	0f af 45 ac          	imul   -0x54(%rbp),%eax
   12988:	48 63 d0             	movslq %eax,%rdx
   1298b:	8b 45 98             	mov    -0x68(%rbp),%eax
   1298e:	48 98                	cltq   
   12990:	48 01 d0             	add    %rdx,%rax
   12993:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   1299a:	00 
   1299b:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1299f:	48 01 d0             	add    %rdx,%rax
   129a2:	66 0f ef c0          	pxor   %xmm0,%xmm0
   129a6:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   129aa:	83 45 98 01          	addl   $0x1,-0x68(%rbp)
   129ae:	8b 45 98             	mov    -0x68(%rbp),%eax
   129b1:	3b 45 ac             	cmp    -0x54(%rbp),%eax
   129b4:	7c cb                	jl     12981 <_Z7paddingPP10particle_tPii+0x2c1>
   129b6:	83 45 94 01          	addl   $0x1,-0x6c(%rbp)
   129ba:	8b 45 94             	mov    -0x6c(%rbp),%eax
   129bd:	3b 45 a4             	cmp    -0x5c(%rbp),%eax
   129c0:	7c b6                	jl     12978 <_Z7paddingPP10particle_tPii+0x2b8>
   129c2:	c7 45 9c 00 00 00 00 	movl   $0x0,-0x64(%rbp)
   129c9:	e9 9b 01 00 00       	jmp    12b69 <_Z7paddingPP10particle_tPii+0x4a9>
   129ce:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   129d2:	48 3b 45 b8          	cmp    -0x48(%rbp),%rax
   129d6:	75 08                	jne    129e0 <_Z7paddingPP10particle_tPii+0x320>
   129d8:	8b 45 9c             	mov    -0x64(%rbp),%eax
   129db:	3b 45 a0             	cmp    -0x60(%rbp),%eax
   129de:	7c 2f                	jl     12a0f <_Z7paddingPP10particle_tPii+0x34f>
   129e0:	8b 45 9c             	mov    -0x64(%rbp),%eax
   129e3:	0f af 45 a8          	imul   -0x58(%rbp),%eax
   129e7:	48 98                	cltq   
   129e9:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   129f0:	00 
   129f1:	8b 45 9c             	mov    -0x64(%rbp),%eax
   129f4:	48 98                	cltq   
   129f6:	48 c1 e0 06          	shl    $0x6,%rax
   129fa:	48 89 c2             	mov    %rax,%rdx
   129fd:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12a01:	48 01 d0             	add    %rdx,%rax
   12a04:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
   12a08:	48 01 ca             	add    %rcx,%rdx
   12a0b:	48 89 50 08          	mov    %rdx,0x8(%rax)
   12a0f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12a13:	48 3b 45 c0          	cmp    -0x40(%rbp),%rax
   12a17:	75 08                	jne    12a21 <_Z7paddingPP10particle_tPii+0x361>
   12a19:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12a1c:	3b 45 a0             	cmp    -0x60(%rbp),%eax
   12a1f:	7c 2f                	jl     12a50 <_Z7paddingPP10particle_tPii+0x390>
   12a21:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12a24:	0f af 45 a8          	imul   -0x58(%rbp),%eax
   12a28:	48 98                	cltq   
   12a2a:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   12a31:	00 
   12a32:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12a35:	48 98                	cltq   
   12a37:	48 c1 e0 06          	shl    $0x6,%rax
   12a3b:	48 89 c2             	mov    %rax,%rdx
   12a3e:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12a42:	48 01 d0             	add    %rdx,%rax
   12a45:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   12a49:	48 01 ca             	add    %rcx,%rdx
   12a4c:	48 89 50 10          	mov    %rdx,0x10(%rax)
   12a50:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   12a54:	48 3b 45 c8          	cmp    -0x38(%rbp),%rax
   12a58:	75 08                	jne    12a62 <_Z7paddingPP10particle_tPii+0x3a2>
   12a5a:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12a5d:	3b 45 a0             	cmp    -0x60(%rbp),%eax
   12a60:	7c 2f                	jl     12a91 <_Z7paddingPP10particle_tPii+0x3d1>
   12a62:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12a65:	0f af 45 a8          	imul   -0x58(%rbp),%eax
   12a69:	48 98                	cltq   
   12a6b:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   12a72:	00 
   12a73:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12a76:	48 98                	cltq   
   12a78:	48 c1 e0 06          	shl    $0x6,%rax
   12a7c:	48 89 c2             	mov    %rax,%rdx
   12a7f:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12a83:	48 01 d0             	add    %rdx,%rax
   12a86:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   12a8a:	48 01 ca             	add    %rcx,%rdx
   12a8d:	48 89 50 18          	mov    %rdx,0x18(%rax)
   12a91:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12a95:	48 3b 45 d0          	cmp    -0x30(%rbp),%rax
   12a99:	75 08                	jne    12aa3 <_Z7paddingPP10particle_tPii+0x3e3>
   12a9b:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12a9e:	3b 45 a0             	cmp    -0x60(%rbp),%eax
   12aa1:	7c 2f                	jl     12ad2 <_Z7paddingPP10particle_tPii+0x412>
   12aa3:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12aa6:	0f af 45 ac          	imul   -0x54(%rbp),%eax
   12aaa:	48 98                	cltq   
   12aac:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   12ab3:	00 
   12ab4:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12ab7:	48 98                	cltq   
   12ab9:	48 c1 e0 06          	shl    $0x6,%rax
   12abd:	48 89 c2             	mov    %rax,%rdx
   12ac0:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12ac4:	48 01 d0             	add    %rdx,%rax
   12ac7:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   12acb:	48 01 ca             	add    %rcx,%rdx
   12ace:	48 89 50 20          	mov    %rdx,0x20(%rax)
   12ad2:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12ad5:	3b 45 a0             	cmp    -0x60(%rbp),%eax
   12ad8:	0f 8c 87 00 00 00    	jl     12b65 <_Z7paddingPP10particle_tPii+0x4a5>
   12ade:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12ae1:	48 98                	cltq   
   12ae3:	48 c1 e0 06          	shl    $0x6,%rax
   12ae7:	48 89 c2             	mov    %rax,%rdx
   12aea:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12aee:	48 01 d0             	add    %rdx,%rax
   12af1:	c6 40 38 00          	movb   $0x0,0x38(%rax)
   12af5:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12af8:	48 63 d0             	movslq %eax,%rdx
   12afb:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12afe:	48 98                	cltq   
   12b00:	48 c1 e0 06          	shl    $0x6,%rax
   12b04:	48 89 c1             	mov    %rax,%rcx
   12b07:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12b0b:	48 01 c8             	add    %rcx,%rax
   12b0e:	48 f7 d2             	not    %rdx
   12b11:	48 89 10             	mov    %rdx,(%rax)
   12b14:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12b17:	48 98                	cltq   
   12b19:	48 c1 e0 06          	shl    $0x6,%rax
   12b1d:	48 89 c2             	mov    %rax,%rdx
   12b20:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12b24:	48 01 c2             	add    %rax,%rdx
   12b27:	8b 45 a8             	mov    -0x58(%rbp),%eax
   12b2a:	89 42 30             	mov    %eax,0x30(%rdx)
   12b2d:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12b30:	48 98                	cltq   
   12b32:	48 c1 e0 06          	shl    $0x6,%rax
   12b36:	48 89 c2             	mov    %rax,%rdx
   12b39:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12b3d:	48 01 c2             	add    %rax,%rdx
   12b40:	8b 45 ac             	mov    -0x54(%rbp),%eax
   12b43:	89 42 34             	mov    %eax,0x34(%rdx)
   12b46:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12b49:	48 98                	cltq   
   12b4b:	48 c1 e0 06          	shl    $0x6,%rax
   12b4f:	48 89 c2             	mov    %rax,%rdx
   12b52:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12b56:	48 01 c2             	add    %rax,%rdx
   12b59:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12b5d:	48 8b 40 28          	mov    0x28(%rax),%rax
   12b61:	48 89 42 28          	mov    %rax,0x28(%rdx)
   12b65:	83 45 9c 01          	addl   $0x1,-0x64(%rbp)
   12b69:	8b 45 9c             	mov    -0x64(%rbp),%eax
   12b6c:	3b 45 a4             	cmp    -0x5c(%rbp),%eax
   12b6f:	0f 8c 59 fe ff ff    	jl     129ce <_Z7paddingPP10particle_tPii+0x30e>
   12b75:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   12b7c:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
   12b80:	48 89 10             	mov    %rdx,(%rax)
   12b83:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   12b8a:	8b 55 a4             	mov    -0x5c(%rbp),%edx
   12b8d:	89 10                	mov    %edx,(%rax)
   12b8f:	90                   	nop
   12b90:	c9                   	leave  
   12b91:	c3                   	ret    

0000000000012b92 <_Z7paddingP17chunk_particles_ti>:
   12b92:	f3 0f 1e fa          	endbr64 
   12b96:	55                   	push   %rbp
   12b97:	48 89 e5             	mov    %rsp,%rbp
   12b9a:	48 83 ec 30          	sub    $0x30,%rsp
   12b9e:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   12ba2:	89 75 d4             	mov    %esi,-0x2c(%rbp)
   12ba5:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   12bac:	00 00 
   12bae:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   12bb2:	31 c0                	xor    %eax,%eax
   12bb4:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12bb8:	48 8b 00             	mov    (%rax),%rax
   12bbb:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   12bbf:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12bc3:	8b 40 08             	mov    0x8(%rax),%eax
   12bc6:	89 45 ec             	mov    %eax,-0x14(%rbp)
   12bc9:	8b 55 d4             	mov    -0x2c(%rbp),%edx
   12bcc:	48 8d 4d ec          	lea    -0x14(%rbp),%rcx
   12bd0:	48 8d 45 f0          	lea    -0x10(%rbp),%rax
   12bd4:	48 89 ce             	mov    %rcx,%rsi
   12bd7:	48 89 c7             	mov    %rax,%rdi
   12bda:	e8 e1 fa ff ff       	call   126c0 <_Z7paddingPP10particle_tPii>
   12bdf:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   12be3:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12be7:	48 89 10             	mov    %rdx,(%rax)
   12bea:	8b 55 ec             	mov    -0x14(%rbp),%edx
   12bed:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12bf1:	89 50 08             	mov    %edx,0x8(%rax)
   12bf4:	90                   	nop
   12bf5:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12bf9:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   12c00:	00 00 
   12c02:	74 05                	je     12c09 <_Z7paddingP17chunk_particles_ti+0x77>
   12c04:	e8 c7 78 ff ff       	call   a4d0 <__stack_chk_fail@plt>
   12c09:	c9                   	leave  
   12c0a:	c3                   	ret    

0000000000012c0b <_Z6shrinkPP10particle_tPi>:
   12c0b:	f3 0f 1e fa          	endbr64 
   12c0f:	55                   	push   %rbp
   12c10:	48 89 e5             	mov    %rsp,%rbp
   12c13:	48 81 ec 90 00 00 00 	sub    $0x90,%rsp
   12c1a:	48 89 bd 78 ff ff ff 	mov    %rdi,-0x88(%rbp)
   12c21:	48 89 b5 70 ff ff ff 	mov    %rsi,-0x90(%rbp)
   12c28:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   12c2f:	48 8b 00             	mov    (%rax),%rax
   12c32:	48 89 45 a8          	mov    %rax,-0x58(%rbp)
   12c36:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   12c3d:	8b 00                	mov    (%rax),%eax
   12c3f:	89 45 94             	mov    %eax,-0x6c(%rbp)
   12c42:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12c46:	48 89 45 98          	mov    %rax,-0x68(%rbp)
   12c4a:	c7 45 80 00 00 00 00 	movl   $0x0,-0x80(%rbp)
   12c51:	eb 18                	jmp    12c6b <_Z6shrinkPP10particle_tPi+0x60>
   12c53:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   12c57:	0f b6 40 38          	movzbl 0x38(%rax),%eax
   12c5b:	83 f0 01             	xor    $0x1,%eax
   12c5e:	84 c0                	test   %al,%al
   12c60:	75 13                	jne    12c75 <_Z6shrinkPP10particle_tPi+0x6a>
   12c62:	48 83 45 98 40       	addq   $0x40,-0x68(%rbp)
   12c67:	83 45 80 01          	addl   $0x1,-0x80(%rbp)
   12c6b:	8b 45 80             	mov    -0x80(%rbp),%eax
   12c6e:	3b 45 94             	cmp    -0x6c(%rbp),%eax
   12c71:	7c e0                	jl     12c53 <_Z6shrinkPP10particle_tPi+0x48>
   12c73:	eb 01                	jmp    12c76 <_Z6shrinkPP10particle_tPi+0x6b>
   12c75:	90                   	nop
   12c76:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12c7a:	48 8b 40 08          	mov    0x8(%rax),%rax
   12c7e:	48 89 45 b0          	mov    %rax,-0x50(%rbp)
   12c82:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12c86:	8b 40 30             	mov    0x30(%rax),%eax
   12c89:	0f af 45 80          	imul   -0x80(%rbp),%eax
   12c8d:	48 98                	cltq   
   12c8f:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12c96:	00 
   12c97:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   12c9b:	48 89 d6             	mov    %rdx,%rsi
   12c9e:	48 89 c7             	mov    %rax,%rdi
   12ca1:	e8 fa 7b ff ff       	call   a8a0 <realloc@plt>
   12ca6:	48 89 45 b8          	mov    %rax,-0x48(%rbp)
   12caa:	48 83 7d b8 00       	cmpq   $0x0,-0x48(%rbp)
   12caf:	75 28                	jne    12cd9 <_Z6shrinkPP10particle_tPi+0xce>
   12cb1:	48 8d 05 b0 40 00 00 	lea    0x40b0(%rip),%rax        # 16d68 <_ZN3MPIL11LOCK_SHAREDE+0x120>
   12cb8:	48 89 c1             	mov    %rax,%rcx
   12cbb:	ba 7e 00 00 00       	mov    $0x7e,%edx
   12cc0:	48 8d 05 b9 3f 00 00 	lea    0x3fb9(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12cc7:	48 89 c6             	mov    %rax,%rsi
   12cca:	48 8d 05 5d 40 00 00 	lea    0x405d(%rip),%rax        # 16d2e <_ZN3MPIL11LOCK_SHAREDE+0xe6>
   12cd1:	48 89 c7             	mov    %rax,%rdi
   12cd4:	e8 57 74 ff ff       	call   a130 <__assert_fail@plt>
   12cd9:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12cdd:	48 8b 40 10          	mov    0x10(%rax),%rax
   12ce1:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   12ce5:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12ce9:	8b 40 30             	mov    0x30(%rax),%eax
   12cec:	0f af 45 80          	imul   -0x80(%rbp),%eax
   12cf0:	48 98                	cltq   
   12cf2:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12cf9:	00 
   12cfa:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   12cfe:	48 89 d6             	mov    %rdx,%rsi
   12d01:	48 89 c7             	mov    %rax,%rdi
   12d04:	e8 97 7b ff ff       	call   a8a0 <realloc@plt>
   12d09:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   12d0d:	48 83 7d c8 00       	cmpq   $0x0,-0x38(%rbp)
   12d12:	75 28                	jne    12d3c <_Z6shrinkPP10particle_tPi+0x131>
   12d14:	48 8d 05 4d 40 00 00 	lea    0x404d(%rip),%rax        # 16d68 <_ZN3MPIL11LOCK_SHAREDE+0x120>
   12d1b:	48 89 c1             	mov    %rax,%rcx
   12d1e:	ba 82 00 00 00       	mov    $0x82,%edx
   12d23:	48 8d 05 56 3f 00 00 	lea    0x3f56(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12d2a:	48 89 c6             	mov    %rax,%rsi
   12d2d:	48 8d 05 0c 40 00 00 	lea    0x400c(%rip),%rax        # 16d40 <_ZN3MPIL11LOCK_SHAREDE+0xf8>
   12d34:	48 89 c7             	mov    %rax,%rdi
   12d37:	e8 f4 73 ff ff       	call   a130 <__assert_fail@plt>
   12d3c:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12d40:	48 8b 40 18          	mov    0x18(%rax),%rax
   12d44:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   12d48:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12d4c:	8b 40 30             	mov    0x30(%rax),%eax
   12d4f:	0f af 45 80          	imul   -0x80(%rbp),%eax
   12d53:	48 98                	cltq   
   12d55:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12d5c:	00 
   12d5d:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   12d61:	48 89 d6             	mov    %rdx,%rsi
   12d64:	48 89 c7             	mov    %rax,%rdi
   12d67:	e8 34 7b ff ff       	call   a8a0 <realloc@plt>
   12d6c:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   12d70:	48 83 7d d8 00       	cmpq   $0x0,-0x28(%rbp)
   12d75:	75 28                	jne    12d9f <_Z6shrinkPP10particle_tPi+0x194>
   12d77:	48 8d 05 ea 3f 00 00 	lea    0x3fea(%rip),%rax        # 16d68 <_ZN3MPIL11LOCK_SHAREDE+0x120>
   12d7e:	48 89 c1             	mov    %rax,%rcx
   12d81:	ba 86 00 00 00       	mov    $0x86,%edx
   12d86:	48 8d 05 f3 3e 00 00 	lea    0x3ef3(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12d8d:	48 89 c6             	mov    %rax,%rsi
   12d90:	48 8d 05 bb 3f 00 00 	lea    0x3fbb(%rip),%rax        # 16d52 <_ZN3MPIL11LOCK_SHAREDE+0x10a>
   12d97:	48 89 c7             	mov    %rax,%rdi
   12d9a:	e8 91 73 ff ff       	call   a130 <__assert_fail@plt>
   12d9f:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12da3:	48 8b 40 20          	mov    0x20(%rax),%rax
   12da7:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   12dab:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12daf:	8b 40 34             	mov    0x34(%rax),%eax
   12db2:	0f af 45 80          	imul   -0x80(%rbp),%eax
   12db6:	48 98                	cltq   
   12db8:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   12dbf:	00 
   12dc0:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   12dc4:	48 89 d6             	mov    %rdx,%rsi
   12dc7:	48 89 c7             	mov    %rax,%rdi
   12dca:	e8 d1 7a ff ff       	call   a8a0 <realloc@plt>
   12dcf:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   12dd3:	48 83 7d e8 00       	cmpq   $0x0,-0x18(%rbp)
   12dd8:	75 28                	jne    12e02 <_Z6shrinkPP10particle_tPi+0x1f7>
   12dda:	48 8d 05 87 3f 00 00 	lea    0x3f87(%rip),%rax        # 16d68 <_ZN3MPIL11LOCK_SHAREDE+0x120>
   12de1:	48 89 c1             	mov    %rax,%rcx
   12de4:	ba 8a 00 00 00       	mov    $0x8a,%edx
   12de9:	48 8d 05 90 3e 00 00 	lea    0x3e90(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12df0:	48 89 c6             	mov    %rax,%rsi
   12df3:	48 8d 05 8e 3f 00 00 	lea    0x3f8e(%rip),%rax        # 16d88 <_ZN3MPIL11LOCK_SHAREDE+0x140>
   12dfa:	48 89 c7             	mov    %rax,%rdi
   12dfd:	e8 2e 73 ff ff       	call   a130 <__assert_fail@plt>
   12e02:	48 8b 45 a8          	mov    -0x58(%rbp),%rax
   12e06:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   12e0a:	8b 45 80             	mov    -0x80(%rbp),%eax
   12e0d:	48 98                	cltq   
   12e0f:	48 c1 e0 06          	shl    $0x6,%rax
   12e13:	48 89 c2             	mov    %rax,%rdx
   12e16:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   12e1a:	48 89 d6             	mov    %rdx,%rsi
   12e1d:	48 89 c7             	mov    %rax,%rdi
   12e20:	e8 7b 7a ff ff       	call   a8a0 <realloc@plt>
   12e25:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   12e29:	48 83 7d f8 00       	cmpq   $0x0,-0x8(%rbp)
   12e2e:	75 28                	jne    12e58 <_Z6shrinkPP10particle_tPi+0x24d>
   12e30:	48 8d 05 31 3f 00 00 	lea    0x3f31(%rip),%rax        # 16d68 <_ZN3MPIL11LOCK_SHAREDE+0x120>
   12e37:	48 89 c1             	mov    %rax,%rcx
   12e3a:	ba 8e 00 00 00       	mov    $0x8e,%edx
   12e3f:	48 8d 05 3a 3e 00 00 	lea    0x3e3a(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   12e46:	48 89 c6             	mov    %rax,%rsi
   12e49:	48 8d 05 ce 3e 00 00 	lea    0x3ece(%rip),%rax        # 16d1e <_ZN3MPIL11LOCK_SHAREDE+0xd6>
   12e50:	48 89 c7             	mov    %rax,%rdi
   12e53:	e8 d8 72 ff ff       	call   a130 <__assert_fail@plt>
   12e58:	48 8b 45 b0          	mov    -0x50(%rbp),%rax
   12e5c:	48 3b 45 b8          	cmp    -0x48(%rbp),%rax
   12e60:	74 47                	je     12ea9 <_Z6shrinkPP10particle_tPi+0x29e>
   12e62:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12e66:	48 89 45 98          	mov    %rax,-0x68(%rbp)
   12e6a:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   12e6e:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
   12e72:	c7 45 84 00 00 00 00 	movl   $0x0,-0x7c(%rbp)
   12e79:	eb 26                	jmp    12ea1 <_Z6shrinkPP10particle_tPi+0x296>
   12e7b:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   12e7f:	48 8b 55 a0          	mov    -0x60(%rbp),%rdx
   12e83:	48 89 50 08          	mov    %rdx,0x8(%rax)
   12e87:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   12e8b:	8b 40 30             	mov    0x30(%rax),%eax
   12e8e:	48 98                	cltq   
   12e90:	48 c1 e0 03          	shl    $0x3,%rax
   12e94:	48 01 45 a0          	add    %rax,-0x60(%rbp)
   12e98:	48 83 45 98 40       	addq   $0x40,-0x68(%rbp)
   12e9d:	83 45 84 01          	addl   $0x1,-0x7c(%rbp)
   12ea1:	8b 45 84             	mov    -0x7c(%rbp),%eax
   12ea4:	3b 45 80             	cmp    -0x80(%rbp),%eax
   12ea7:	7c d2                	jl     12e7b <_Z6shrinkPP10particle_tPi+0x270>
   12ea9:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   12ead:	48 3b 45 c8          	cmp    -0x38(%rbp),%rax
   12eb1:	74 47                	je     12efa <_Z6shrinkPP10particle_tPi+0x2ef>
   12eb3:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12eb7:	48 89 45 98          	mov    %rax,-0x68(%rbp)
   12ebb:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   12ebf:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
   12ec3:	c7 45 88 00 00 00 00 	movl   $0x0,-0x78(%rbp)
   12eca:	eb 26                	jmp    12ef2 <_Z6shrinkPP10particle_tPi+0x2e7>
   12ecc:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   12ed0:	48 8b 55 a0          	mov    -0x60(%rbp),%rdx
   12ed4:	48 89 50 10          	mov    %rdx,0x10(%rax)
   12ed8:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   12edc:	8b 40 30             	mov    0x30(%rax),%eax
   12edf:	48 98                	cltq   
   12ee1:	48 c1 e0 03          	shl    $0x3,%rax
   12ee5:	48 01 45 a0          	add    %rax,-0x60(%rbp)
   12ee9:	48 83 45 98 40       	addq   $0x40,-0x68(%rbp)
   12eee:	83 45 88 01          	addl   $0x1,-0x78(%rbp)
   12ef2:	8b 45 88             	mov    -0x78(%rbp),%eax
   12ef5:	3b 45 80             	cmp    -0x80(%rbp),%eax
   12ef8:	7c d2                	jl     12ecc <_Z6shrinkPP10particle_tPi+0x2c1>
   12efa:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   12efe:	48 3b 45 d8          	cmp    -0x28(%rbp),%rax
   12f02:	74 47                	je     12f4b <_Z6shrinkPP10particle_tPi+0x340>
   12f04:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12f08:	48 89 45 98          	mov    %rax,-0x68(%rbp)
   12f0c:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12f10:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
   12f14:	c7 45 8c 00 00 00 00 	movl   $0x0,-0x74(%rbp)
   12f1b:	eb 26                	jmp    12f43 <_Z6shrinkPP10particle_tPi+0x338>
   12f1d:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   12f21:	48 8b 55 a0          	mov    -0x60(%rbp),%rdx
   12f25:	48 89 50 18          	mov    %rdx,0x18(%rax)
   12f29:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   12f2d:	8b 40 30             	mov    0x30(%rax),%eax
   12f30:	48 98                	cltq   
   12f32:	48 c1 e0 03          	shl    $0x3,%rax
   12f36:	48 01 45 a0          	add    %rax,-0x60(%rbp)
   12f3a:	48 83 45 98 40       	addq   $0x40,-0x68(%rbp)
   12f3f:	83 45 8c 01          	addl   $0x1,-0x74(%rbp)
   12f43:	8b 45 8c             	mov    -0x74(%rbp),%eax
   12f46:	3b 45 80             	cmp    -0x80(%rbp),%eax
   12f49:	7c d2                	jl     12f1d <_Z6shrinkPP10particle_tPi+0x312>
   12f4b:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   12f4f:	48 3b 45 e8          	cmp    -0x18(%rbp),%rax
   12f53:	74 47                	je     12f9c <_Z6shrinkPP10particle_tPi+0x391>
   12f55:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   12f59:	48 89 45 98          	mov    %rax,-0x68(%rbp)
   12f5d:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   12f61:	48 89 45 a0          	mov    %rax,-0x60(%rbp)
   12f65:	c7 45 90 00 00 00 00 	movl   $0x0,-0x70(%rbp)
   12f6c:	eb 26                	jmp    12f94 <_Z6shrinkPP10particle_tPi+0x389>
   12f6e:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   12f72:	48 8b 55 a0          	mov    -0x60(%rbp),%rdx
   12f76:	48 89 50 20          	mov    %rdx,0x20(%rax)
   12f7a:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   12f7e:	8b 40 34             	mov    0x34(%rax),%eax
   12f81:	48 98                	cltq   
   12f83:	48 c1 e0 03          	shl    $0x3,%rax
   12f87:	48 01 45 a0          	add    %rax,-0x60(%rbp)
   12f8b:	48 83 45 98 40       	addq   $0x40,-0x68(%rbp)
   12f90:	83 45 90 01          	addl   $0x1,-0x70(%rbp)
   12f94:	8b 45 90             	mov    -0x70(%rbp),%eax
   12f97:	3b 45 80             	cmp    -0x80(%rbp),%eax
   12f9a:	7c d2                	jl     12f6e <_Z6shrinkPP10particle_tPi+0x363>
   12f9c:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   12fa3:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   12fa7:	48 89 10             	mov    %rdx,(%rax)
   12faa:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   12fb1:	8b 55 80             	mov    -0x80(%rbp),%edx
   12fb4:	89 10                	mov    %edx,(%rax)
   12fb6:	90                   	nop
   12fb7:	c9                   	leave  
   12fb8:	c3                   	ret    

0000000000012fb9 <_Z6shrinkP17chunk_particles_t>:
   12fb9:	f3 0f 1e fa          	endbr64 
   12fbd:	55                   	push   %rbp
   12fbe:	48 89 e5             	mov    %rsp,%rbp
   12fc1:	48 83 ec 30          	sub    $0x30,%rsp
   12fc5:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   12fc9:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   12fd0:	00 00 
   12fd2:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   12fd6:	31 c0                	xor    %eax,%eax
   12fd8:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12fdc:	48 8b 00             	mov    (%rax),%rax
   12fdf:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   12fe3:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   12fe7:	8b 40 08             	mov    0x8(%rax),%eax
   12fea:	89 45 ec             	mov    %eax,-0x14(%rbp)
   12fed:	48 8d 55 ec          	lea    -0x14(%rbp),%rdx
   12ff1:	48 8d 45 f0          	lea    -0x10(%rbp),%rax
   12ff5:	48 89 d6             	mov    %rdx,%rsi
   12ff8:	48 89 c7             	mov    %rax,%rdi
   12ffb:	e8 0b fc ff ff       	call   12c0b <_Z6shrinkPP10particle_tPi>
   13000:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   13004:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   13008:	48 89 10             	mov    %rdx,(%rax)
   1300b:	8b 55 ec             	mov    -0x14(%rbp),%edx
   1300e:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   13012:	89 50 08             	mov    %edx,0x8(%rax)
   13015:	90                   	nop
   13016:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1301a:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   13021:	00 00 
   13023:	74 05                	je     1302a <_Z6shrinkP17chunk_particles_t+0x71>
   13025:	e8 a6 74 ff ff       	call   a4d0 <__stack_chk_fail@plt>
   1302a:	c9                   	leave  
   1302b:	c3                   	ret    

000000000001302c <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t>:
   1302c:	f3 0f 1e fa          	endbr64 
   13030:	55                   	push   %rbp
   13031:	48 89 e5             	mov    %rsp,%rbp
   13034:	48 81 ec a0 00 00 00 	sub    $0xa0,%rsp
   1303b:	48 89 7d 88          	mov    %rdi,-0x78(%rbp)
   1303f:	48 89 75 80          	mov    %rsi,-0x80(%rbp)
   13043:	48 89 95 78 ff ff ff 	mov    %rdx,-0x88(%rbp)
   1304a:	48 89 8d 70 ff ff ff 	mov    %rcx,-0x90(%rbp)
   13051:	4c 89 85 68 ff ff ff 	mov    %r8,-0x98(%rbp)
   13058:	44 89 8d 64 ff ff ff 	mov    %r9d,-0x9c(%rbp)
   1305f:	c6 45 9f 00          	movb   $0x0,-0x61(%rbp)
   13063:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   1306a:	48 85 c0             	test   %rax,%rax
   1306d:	75 22                	jne    13091 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x65>
   1306f:	c6 45 9f 01          	movb   $0x1,-0x61(%rbp)
   13073:	48 8d 85 68 ff ff ff 	lea    -0x98(%rbp),%rax
   1307a:	b9 01 00 00 00       	mov    $0x1,%ecx
   1307f:	ba 01 00 00 00       	mov    $0x1,%edx
   13084:	be 01 00 00 00       	mov    $0x1,%esi
   13089:	48 89 c7             	mov    %rax,%rdi
   1308c:	e8 85 ed ff ff       	call   11e16 <_Z15alloc_particlesPP10particle_tiii>
   13091:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13095:	48 89 45 b8          	mov    %rax,-0x48(%rbp)
   13099:	48 8b 45 10          	mov    0x10(%rbp),%rax
   1309d:	48 83 c0 04          	add    $0x4,%rax
   130a1:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   130a5:	48 8b 45 10          	mov    0x10(%rbp),%rax
   130a9:	48 83 c0 08          	add    $0x8,%rax
   130ad:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   130b1:	48 8b 45 10          	mov    0x10(%rbp),%rax
   130b5:	48 83 c0 0c          	add    $0xc,%rax
   130b9:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   130bd:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   130c1:	8b 00                	mov    (%rax),%eax
   130c3:	85 c0                	test   %eax,%eax
   130c5:	0f 85 84 00 00 00    	jne    1314f <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x123>
   130cb:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   130cf:	8b 00                	mov    (%rax),%eax
   130d1:	85 c0                	test   %eax,%eax
   130d3:	75 7a                	jne    1314f <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x123>
   130d5:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   130d9:	8b 08                	mov    (%rax),%ecx
   130db:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   130e1:	99                   	cltd   
   130e2:	f7 f9                	idiv   %ecx
   130e4:	89 d0                	mov    %edx,%eax
   130e6:	85 c0                	test   %eax,%eax
   130e8:	74 28                	je     13112 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0xe6>
   130ea:	48 8d 05 af 3c 00 00 	lea    0x3caf(%rip),%rax        # 16da0 <_ZN3MPIL11LOCK_SHAREDE+0x158>
   130f1:	48 89 c1             	mov    %rax,%rcx
   130f4:	ba d3 00 00 00       	mov    $0xd3,%edx
   130f9:	48 8d 05 80 3b 00 00 	lea    0x3b80(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   13100:	48 89 c6             	mov    %rax,%rsi
   13103:	48 8d 05 ea 3c 00 00 	lea    0x3cea(%rip),%rax        # 16df4 <_ZN3MPIL11LOCK_SHAREDE+0x1ac>
   1310a:	48 89 c7             	mov    %rax,%rdi
   1310d:	e8 1e 70 ff ff       	call   a130 <__assert_fail@plt>
   13112:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   13116:	8b 08                	mov    (%rax),%ecx
   13118:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   1311e:	99                   	cltd   
   1311f:	f7 f9                	idiv   %ecx
   13121:	89 d0                	mov    %edx,%eax
   13123:	85 c0                	test   %eax,%eax
   13125:	74 28                	je     1314f <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x123>
   13127:	48 8d 05 72 3c 00 00 	lea    0x3c72(%rip),%rax        # 16da0 <_ZN3MPIL11LOCK_SHAREDE+0x158>
   1312e:	48 89 c1             	mov    %rax,%rcx
   13131:	ba d4 00 00 00       	mov    $0xd4,%edx
   13136:	48 8d 05 43 3b 00 00 	lea    0x3b43(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   1313d:	48 89 c6             	mov    %rax,%rsi
   13140:	48 8d 05 c8 3c 00 00 	lea    0x3cc8(%rip),%rax        # 16e0f <_ZN3MPIL11LOCK_SHAREDE+0x1c7>
   13147:	48 89 c7             	mov    %rax,%rdi
   1314a:	e8 e1 6f ff ff       	call   a130 <__assert_fail@plt>
   1314f:	bf 10 00 00 00       	mov    $0x10,%edi
   13154:	e8 c7 74 ff ff       	call   a620 <malloc@plt>
   13159:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   1315d:	48 83 7d d8 00       	cmpq   $0x0,-0x28(%rbp)
   13162:	75 28                	jne    1318c <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x160>
   13164:	48 8d 05 35 3c 00 00 	lea    0x3c35(%rip),%rax        # 16da0 <_ZN3MPIL11LOCK_SHAREDE+0x158>
   1316b:	48 89 c1             	mov    %rax,%rcx
   1316e:	ba d8 00 00 00       	mov    $0xd8,%edx
   13173:	48 8d 05 06 3b 00 00 	lea    0x3b06(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   1317a:	48 89 c6             	mov    %rax,%rsi
   1317d:	48 8d 05 ac 3c 00 00 	lea    0x3cac(%rip),%rax        # 16e30 <_ZN3MPIL11LOCK_SHAREDE+0x1e8>
   13184:	48 89 c7             	mov    %rax,%rdi
   13187:	e8 a4 6f ff ff       	call   a130 <__assert_fail@plt>
   1318c:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   13190:	8b 00                	mov    (%rax),%eax
   13192:	85 c0                	test   %eax,%eax
   13194:	75 62                	jne    131f8 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x1cc>
   13196:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   1319a:	8b 00                	mov    (%rax),%eax
   1319c:	85 c0                	test   %eax,%eax
   1319e:	75 58                	jne    131f8 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x1cc>
   131a0:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   131a4:	8b 38                	mov    (%rax),%edi
   131a6:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   131ac:	99                   	cltd   
   131ad:	f7 ff                	idiv   %edi
   131af:	89 c2                	mov    %eax,%edx
   131b1:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   131b5:	89 10                	mov    %edx,(%rax)
   131b7:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   131bb:	8b 38                	mov    (%rax),%edi
   131bd:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   131c1:	48 8d 48 04          	lea    0x4(%rax),%rcx
   131c5:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   131cb:	99                   	cltd   
   131cc:	f7 ff                	idiv   %edi
   131ce:	89 01                	mov    %eax,(%rcx)
   131d0:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   131d7:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
   131db:	48 83 c2 08          	add    $0x8,%rdx
   131df:	8b 40 30             	mov    0x30(%rax),%eax
   131e2:	89 02                	mov    %eax,(%rdx)
   131e4:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   131eb:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
   131ef:	48 83 c2 0c          	add    $0xc,%rdx
   131f3:	8b 40 34             	mov    0x34(%rax),%eax
   131f6:	89 02                	mov    %eax,(%rdx)
   131f8:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   131fc:	4c 8d 05 1d a8 00 00 	lea    0xa81d(%rip),%r8        # 1da20 <ompi_mpi_comm_world@@Base>
   13203:	b9 00 00 00 00       	mov    $0x0,%ecx
   13208:	48 8d 15 11 ae 00 00 	lea    0xae11(%rip),%rdx        # 1e020 <ompi_mpi_int@@Base>
   1320f:	be 04 00 00 00       	mov    $0x4,%esi
   13214:	48 89 c7             	mov    %rax,%rdi
   13217:	e8 34 6c ff ff       	call   9e50 <MPI_Bcast@plt>
   1321c:	48 8d 05 fd a7 00 00 	lea    0xa7fd(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   13223:	48 89 c7             	mov    %rax,%rdi
   13226:	e8 95 6f ff ff       	call   a1c0 <MPI_Barrier@plt>
   1322b:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   1322f:	8b 10                	mov    (%rax),%edx
   13231:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   13235:	89 10                	mov    %edx,(%rax)
   13237:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   1323b:	8b 50 04             	mov    0x4(%rax),%edx
   1323e:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   13245:	89 10                	mov    %edx,(%rax)
   13247:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   1324b:	8b 40 08             	mov    0x8(%rax),%eax
   1324e:	89 45 b0             	mov    %eax,-0x50(%rbp)
   13251:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   13255:	8b 40 0c             	mov    0xc(%rax),%eax
   13258:	89 45 b4             	mov    %eax,-0x4c(%rbp)
   1325b:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   1325f:	8b 30                	mov    (%rax),%esi
   13261:	8b 4d b4             	mov    -0x4c(%rbp),%ecx
   13264:	8b 55 b0             	mov    -0x50(%rbp),%edx
   13267:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   1326b:	48 89 c7             	mov    %rax,%rdi
   1326e:	e8 a3 eb ff ff       	call   11e16 <_Z15alloc_particlesPP10particle_tiii>
   13273:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   1327a:	8b 30                	mov    (%rax),%esi
   1327c:	8b 4d b4             	mov    -0x4c(%rbp),%ecx
   1327f:	8b 55 b0             	mov    -0x50(%rbp),%edx
   13282:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13289:	48 89 c7             	mov    %rax,%rdi
   1328c:	e8 85 eb ff ff       	call   11e16 <_Z15alloc_particlesPP10particle_tiii>
   13291:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   13295:	8b 00                	mov    (%rax),%eax
   13297:	85 c0                	test   %eax,%eax
   13299:	0f 85 be 04 00 00    	jne    1375d <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x731>
   1329f:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   132a3:	48 8b 00             	mov    (%rax),%rax
   132a6:	48 8b 40 08          	mov    0x8(%rax),%rax
   132aa:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   132ae:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   132b2:	48 8b 00             	mov    (%rax),%rax
   132b5:	48 8b 40 10          	mov    0x10(%rax),%rax
   132b9:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   132bd:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   132c1:	48 8b 00             	mov    (%rax),%rax
   132c4:	48 8b 40 18          	mov    0x18(%rax),%rax
   132c8:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   132cc:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   132d0:	48 8b 00             	mov    (%rax),%rax
   132d3:	48 8b 40 20          	mov    0x20(%rax),%rax
   132d7:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   132db:	48 8b 4d f0          	mov    -0x10(%rbp),%rcx
   132df:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   132e3:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   132e7:	48 83 ec 08          	sub    $0x8,%rsp
   132eb:	ff 75 f8             	push   -0x8(%rbp)
   132ee:	49 89 c9             	mov    %rcx,%r9
   132f1:	49 89 d0             	mov    %rdx,%r8
   132f4:	48 89 c1             	mov    %rax,%rcx
   132f7:	ba f3 00 00 00       	mov    $0xf3,%edx
   132fc:	48 8d 05 7d 39 00 00 	lea    0x397d(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   13303:	48 89 c6             	mov    %rax,%rsi
   13306:	48 8d 05 45 3b 00 00 	lea    0x3b45(%rip),%rax        # 16e52 <_ZN3MPIL11LOCK_SHAREDE+0x20a>
   1330d:	48 89 c7             	mov    %rax,%rdi
   13310:	b8 00 00 00 00       	mov    $0x0,%eax
   13315:	e8 a6 6a ff ff       	call   9dc0 <printf@plt>
   1331a:	48 83 c4 10          	add    $0x10,%rsp
   1331e:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13322:	48 8b 40 20          	mov    0x20(%rax),%rax
   13326:	48 89 c7             	mov    %rax,%rdi
   13329:	e8 92 6e ff ff       	call   a1c0 <MPI_Barrier@plt>
   1332e:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   13332:	8b 30                	mov    (%rax),%esi
   13334:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13338:	8b 50 04             	mov    0x4(%rax),%edx
   1333b:	48 8b 45 10          	mov    0x10(%rbp),%rax
   1333f:	8b 00                	mov    (%rax),%eax
   13341:	48 83 ec 08          	sub    $0x8,%rsp
   13345:	8b 8d 64 ff ff ff    	mov    -0x9c(%rbp),%ecx
   1334b:	51                   	push   %rcx
   1334c:	41 89 f1             	mov    %esi,%r9d
   1334f:	41 b8 f6 00 00 00    	mov    $0xf6,%r8d
   13355:	48 8d 0d 24 39 00 00 	lea    0x3924(%rip),%rcx        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   1335c:	89 c6                	mov    %eax,%esi
   1335e:	48 8d 05 05 3b 00 00 	lea    0x3b05(%rip),%rax        # 16e6a <_ZN3MPIL11LOCK_SHAREDE+0x222>
   13365:	48 89 c7             	mov    %rax,%rdi
   13368:	b8 00 00 00 00       	mov    $0x0,%eax
   1336d:	e8 4e 6a ff ff       	call   9dc0 <printf@plt>
   13372:	48 83 c4 10          	add    $0x10,%rsp
   13376:	48 8b 45 10          	mov    0x10(%rbp),%rax
   1337a:	48 8b 70 20          	mov    0x20(%rax),%rsi
   1337e:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   13382:	8b 00                	mov    (%rax),%eax
   13384:	48 98                	cltq   
   13386:	c1 e0 06             	shl    $0x6,%eax
   13389:	41 89 c0             	mov    %eax,%r8d
   1338c:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   13390:	48 8b 08             	mov    (%rax),%rcx
   13393:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   13399:	48 98                	cltq   
   1339b:	48 c1 e0 06          	shl    $0x6,%rax
   1339f:	48 89 c7             	mov    %rax,%rdi
   133a2:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   133a6:	8b 00                	mov    (%rax),%eax
   133a8:	4c 63 d8             	movslq %eax,%r11
   133ab:	48 89 f8             	mov    %rdi,%rax
   133ae:	ba 00 00 00 00       	mov    $0x0,%edx
   133b3:	49 f7 f3             	div    %r11
   133b6:	89 c7                	mov    %eax,%edi
   133b8:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   133bf:	56                   	push   %rsi
   133c0:	6a 00                	push   $0x0
   133c2:	4c 8d 0d 57 a4 00 00 	lea    0xa457(%rip),%r9        # 1d820 <ompi_mpi_byte@@Base>
   133c9:	48 8d 15 50 a4 00 00 	lea    0xa450(%rip),%rdx        # 1d820 <ompi_mpi_byte@@Base>
   133d0:	89 fe                	mov    %edi,%esi
   133d2:	48 89 c7             	mov    %rax,%rdi
   133d5:	e8 76 6d ff ff       	call   a150 <MPI_Scatter@plt>
   133da:	48 83 c4 10          	add    $0x10,%rsp
   133de:	c7 45 a0 00 00 00 00 	movl   $0x0,-0x60(%rbp)
   133e5:	e9 c0 00 00 00       	jmp    134aa <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x47e>
   133ea:	8b 45 a0             	mov    -0x60(%rbp),%eax
   133ed:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   133f1:	48 98                	cltq   
   133f3:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   133fa:	00 
   133fb:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   133ff:	48 8b 10             	mov    (%rax),%rdx
   13402:	8b 45 a0             	mov    -0x60(%rbp),%eax
   13405:	48 98                	cltq   
   13407:	48 c1 e0 06          	shl    $0x6,%rax
   1340b:	48 01 d0             	add    %rdx,%rax
   1340e:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
   13412:	48 01 ca             	add    %rcx,%rdx
   13415:	48 89 50 08          	mov    %rdx,0x8(%rax)
   13419:	8b 45 a0             	mov    -0x60(%rbp),%eax
   1341c:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13420:	48 98                	cltq   
   13422:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13429:	00 
   1342a:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   1342e:	48 8b 10             	mov    (%rax),%rdx
   13431:	8b 45 a0             	mov    -0x60(%rbp),%eax
   13434:	48 98                	cltq   
   13436:	48 c1 e0 06          	shl    $0x6,%rax
   1343a:	48 01 d0             	add    %rdx,%rax
   1343d:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   13441:	48 01 ca             	add    %rcx,%rdx
   13444:	48 89 50 10          	mov    %rdx,0x10(%rax)
   13448:	8b 45 a0             	mov    -0x60(%rbp),%eax
   1344b:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   1344f:	48 98                	cltq   
   13451:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13458:	00 
   13459:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   1345d:	48 8b 10             	mov    (%rax),%rdx
   13460:	8b 45 a0             	mov    -0x60(%rbp),%eax
   13463:	48 98                	cltq   
   13465:	48 c1 e0 06          	shl    $0x6,%rax
   13469:	48 01 d0             	add    %rdx,%rax
   1346c:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   13470:	48 01 ca             	add    %rcx,%rdx
   13473:	48 89 50 18          	mov    %rdx,0x18(%rax)
   13477:	8b 45 a0             	mov    -0x60(%rbp),%eax
   1347a:	0f af 45 b4          	imul   -0x4c(%rbp),%eax
   1347e:	48 98                	cltq   
   13480:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13487:	00 
   13488:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   1348c:	48 8b 10             	mov    (%rax),%rdx
   1348f:	8b 45 a0             	mov    -0x60(%rbp),%eax
   13492:	48 98                	cltq   
   13494:	48 c1 e0 06          	shl    $0x6,%rax
   13498:	48 01 d0             	add    %rdx,%rax
   1349b:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   1349f:	48 01 ca             	add    %rcx,%rdx
   134a2:	48 89 50 20          	mov    %rdx,0x20(%rax)
   134a6:	83 45 a0 01          	addl   $0x1,-0x60(%rbp)
   134aa:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   134ae:	8b 00                	mov    (%rax),%eax
   134b0:	39 45 a0             	cmp    %eax,-0x60(%rbp)
   134b3:	0f 8c 31 ff ff ff    	jl     133ea <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x3be>
   134b9:	48 8b 45 10          	mov    0x10(%rbp),%rax
   134bd:	48 8b 40 20          	mov    0x20(%rax),%rax
   134c1:	48 89 c7             	mov    %rax,%rdi
   134c4:	e8 f7 6c ff ff       	call   a1c0 <MPI_Barrier@plt>
   134c9:	48 8b 45 10          	mov    0x10(%rbp),%rax
   134cd:	48 8b 48 20          	mov    0x20(%rax),%rcx
   134d1:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   134d5:	8b 00                	mov    (%rax),%eax
   134d7:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   134db:	89 c6                	mov    %eax,%esi
   134dd:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   134e3:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   134e7:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
   134eb:	8b 3a                	mov    (%rdx),%edi
   134ed:	99                   	cltd   
   134ee:	f7 ff                	idiv   %edi
   134f0:	89 c7                	mov    %eax,%edi
   134f2:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   134f9:	48 8b 40 08          	mov    0x8(%rax),%rax
   134fd:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
   13501:	51                   	push   %rcx
   13502:	6a 00                	push   $0x0
   13504:	4c 8d 0d 15 a7 00 00 	lea    0xa715(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   1350b:	41 89 f0             	mov    %esi,%r8d
   1350e:	48 89 d1             	mov    %rdx,%rcx
   13511:	48 8d 15 08 a7 00 00 	lea    0xa708(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13518:	89 fe                	mov    %edi,%esi
   1351a:	48 89 c7             	mov    %rax,%rdi
   1351d:	e8 2e 6c ff ff       	call   a150 <MPI_Scatter@plt>
   13522:	48 83 c4 10          	add    $0x10,%rsp
   13526:	48 8b 45 10          	mov    0x10(%rbp),%rax
   1352a:	48 8b 40 20          	mov    0x20(%rax),%rax
   1352e:	48 89 c7             	mov    %rax,%rdi
   13531:	e8 8a 6c ff ff       	call   a1c0 <MPI_Barrier@plt>
   13536:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   1353a:	8b 30                	mov    (%rax),%esi
   1353c:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13540:	8b 50 04             	mov    0x4(%rax),%edx
   13543:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13547:	8b 00                	mov    (%rax),%eax
   13549:	48 83 ec 08          	sub    $0x8,%rsp
   1354d:	8b 8d 64 ff ff ff    	mov    -0x9c(%rbp),%ecx
   13553:	51                   	push   %rcx
   13554:	41 89 f1             	mov    %esi,%r9d
   13557:	41 b8 17 01 00 00    	mov    $0x117,%r8d
   1355d:	48 8d 0d 1c 37 00 00 	lea    0x371c(%rip),%rcx        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   13564:	89 c6                	mov    %eax,%esi
   13566:	48 8d 05 fd 38 00 00 	lea    0x38fd(%rip),%rax        # 16e6a <_ZN3MPIL11LOCK_SHAREDE+0x222>
   1356d:	48 89 c7             	mov    %rax,%rdi
   13570:	b8 00 00 00 00       	mov    $0x0,%eax
   13575:	e8 46 68 ff ff       	call   9dc0 <printf@plt>
   1357a:	48 83 c4 10          	add    $0x10,%rsp
   1357e:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13582:	48 8b 70 20          	mov    0x20(%rax),%rsi
   13586:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   1358a:	8b 00                	mov    (%rax),%eax
   1358c:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13590:	89 c7                	mov    %eax,%edi
   13592:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   13596:	48 8b 00             	mov    (%rax),%rax
   13599:	48 8b 48 10          	mov    0x10(%rax),%rcx
   1359d:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   135a3:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   135a7:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
   135ab:	44 8b 1a             	mov    (%rdx),%r11d
   135ae:	99                   	cltd   
   135af:	41 f7 fb             	idiv   %r11d
   135b2:	41 89 c2             	mov    %eax,%r10d
   135b5:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   135bc:	48 8b 40 10          	mov    0x10(%rax),%rax
   135c0:	56                   	push   %rsi
   135c1:	6a 00                	push   $0x0
   135c3:	4c 8d 0d 56 a6 00 00 	lea    0xa656(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   135ca:	41 89 f8             	mov    %edi,%r8d
   135cd:	48 8d 15 4c a6 00 00 	lea    0xa64c(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   135d4:	44 89 d6             	mov    %r10d,%esi
   135d7:	48 89 c7             	mov    %rax,%rdi
   135da:	e8 71 6b ff ff       	call   a150 <MPI_Scatter@plt>
   135df:	48 83 c4 10          	add    $0x10,%rsp
   135e3:	48 8b 45 10          	mov    0x10(%rbp),%rax
   135e7:	48 8b 40 20          	mov    0x20(%rax),%rax
   135eb:	48 89 c7             	mov    %rax,%rdi
   135ee:	e8 cd 6b ff ff       	call   a1c0 <MPI_Barrier@plt>
   135f3:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   135f7:	8b 30                	mov    (%rax),%esi
   135f9:	48 8b 45 10          	mov    0x10(%rbp),%rax
   135fd:	8b 50 04             	mov    0x4(%rax),%edx
   13600:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13604:	8b 00                	mov    (%rax),%eax
   13606:	48 83 ec 08          	sub    $0x8,%rsp
   1360a:	8b 8d 64 ff ff ff    	mov    -0x9c(%rbp),%ecx
   13610:	51                   	push   %rcx
   13611:	41 89 f1             	mov    %esi,%r9d
   13614:	41 b8 24 01 00 00    	mov    $0x124,%r8d
   1361a:	48 8d 0d 5f 36 00 00 	lea    0x365f(%rip),%rcx        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   13621:	89 c6                	mov    %eax,%esi
   13623:	48 8d 05 40 38 00 00 	lea    0x3840(%rip),%rax        # 16e6a <_ZN3MPIL11LOCK_SHAREDE+0x222>
   1362a:	48 89 c7             	mov    %rax,%rdi
   1362d:	b8 00 00 00 00       	mov    $0x0,%eax
   13632:	e8 89 67 ff ff       	call   9dc0 <printf@plt>
   13637:	48 83 c4 10          	add    $0x10,%rsp
   1363b:	48 8b 45 10          	mov    0x10(%rbp),%rax
   1363f:	48 8b 70 20          	mov    0x20(%rax),%rsi
   13643:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   13647:	8b 00                	mov    (%rax),%eax
   13649:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   1364d:	89 c7                	mov    %eax,%edi
   1364f:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   13653:	48 8b 00             	mov    (%rax),%rax
   13656:	48 8b 48 18          	mov    0x18(%rax),%rcx
   1365a:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   13660:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13664:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
   13668:	44 8b 1a             	mov    (%rdx),%r11d
   1366b:	99                   	cltd   
   1366c:	41 f7 fb             	idiv   %r11d
   1366f:	41 89 c2             	mov    %eax,%r10d
   13672:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   13679:	48 8b 40 18          	mov    0x18(%rax),%rax
   1367d:	56                   	push   %rsi
   1367e:	6a 00                	push   $0x0
   13680:	4c 8d 0d 99 a5 00 00 	lea    0xa599(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   13687:	41 89 f8             	mov    %edi,%r8d
   1368a:	48 8d 15 8f a5 00 00 	lea    0xa58f(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13691:	44 89 d6             	mov    %r10d,%esi
   13694:	48 89 c7             	mov    %rax,%rdi
   13697:	e8 b4 6a ff ff       	call   a150 <MPI_Scatter@plt>
   1369c:	48 83 c4 10          	add    $0x10,%rsp
   136a0:	48 8b 45 10          	mov    0x10(%rbp),%rax
   136a4:	48 8b 40 20          	mov    0x20(%rax),%rax
   136a8:	48 89 c7             	mov    %rax,%rdi
   136ab:	e8 10 6b ff ff       	call   a1c0 <MPI_Barrier@plt>
   136b0:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   136b4:	8b 30                	mov    (%rax),%esi
   136b6:	48 8b 45 10          	mov    0x10(%rbp),%rax
   136ba:	8b 50 04             	mov    0x4(%rax),%edx
   136bd:	48 8b 45 10          	mov    0x10(%rbp),%rax
   136c1:	8b 00                	mov    (%rax),%eax
   136c3:	48 83 ec 08          	sub    $0x8,%rsp
   136c7:	8b 8d 64 ff ff ff    	mov    -0x9c(%rbp),%ecx
   136cd:	51                   	push   %rcx
   136ce:	41 89 f1             	mov    %esi,%r9d
   136d1:	41 b8 31 01 00 00    	mov    $0x131,%r8d
   136d7:	48 8d 0d a2 35 00 00 	lea    0x35a2(%rip),%rcx        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   136de:	89 c6                	mov    %eax,%esi
   136e0:	48 8d 05 83 37 00 00 	lea    0x3783(%rip),%rax        # 16e6a <_ZN3MPIL11LOCK_SHAREDE+0x222>
   136e7:	48 89 c7             	mov    %rax,%rdi
   136ea:	b8 00 00 00 00       	mov    $0x0,%eax
   136ef:	e8 cc 66 ff ff       	call   9dc0 <printf@plt>
   136f4:	48 83 c4 10          	add    $0x10,%rsp
   136f8:	48 8b 45 10          	mov    0x10(%rbp),%rax
   136fc:	48 8b 70 20          	mov    0x20(%rax),%rsi
   13700:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   13704:	8b 00                	mov    (%rax),%eax
   13706:	0f af 45 b4          	imul   -0x4c(%rbp),%eax
   1370a:	89 c7                	mov    %eax,%edi
   1370c:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   13710:	48 8b 00             	mov    (%rax),%rax
   13713:	48 8b 48 20          	mov    0x20(%rax),%rcx
   13717:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   1371d:	0f af 45 b4          	imul   -0x4c(%rbp),%eax
   13721:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
   13725:	44 8b 1a             	mov    (%rdx),%r11d
   13728:	99                   	cltd   
   13729:	41 f7 fb             	idiv   %r11d
   1372c:	41 89 c2             	mov    %eax,%r10d
   1372f:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   13736:	48 8b 40 20          	mov    0x20(%rax),%rax
   1373a:	56                   	push   %rsi
   1373b:	6a 00                	push   $0x0
   1373d:	4c 8d 0d dc a4 00 00 	lea    0xa4dc(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   13744:	41 89 f8             	mov    %edi,%r8d
   13747:	48 8d 15 d2 a4 00 00 	lea    0xa4d2(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   1374e:	44 89 d6             	mov    %r10d,%esi
   13751:	48 89 c7             	mov    %rax,%rdi
   13754:	e8 f7 69 ff ff       	call   a150 <MPI_Scatter@plt>
   13759:	48 83 c4 10          	add    $0x10,%rsp
   1375d:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   13761:	48 8b 00             	mov    (%rax),%rax
   13764:	48 8b 40 08          	mov    0x8(%rax),%rax
   13768:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   1376c:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   13770:	48 8b 00             	mov    (%rax),%rax
   13773:	48 8b 40 10          	mov    0x10(%rax),%rax
   13777:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   1377b:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   1377f:	48 8b 00             	mov    (%rax),%rax
   13782:	48 8b 40 18          	mov    0x18(%rax),%rax
   13786:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   1378a:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   1378e:	48 8b 00             	mov    (%rax),%rax
   13791:	48 8b 40 20          	mov    0x20(%rax),%rax
   13795:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   13799:	48 8b 45 10          	mov    0x10(%rbp),%rax
   1379d:	48 8b 50 10          	mov    0x10(%rax),%rdx
   137a1:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   137a5:	8b 00                	mov    (%rax),%eax
   137a7:	48 98                	cltq   
   137a9:	c1 e0 06             	shl    $0x6,%eax
   137ac:	89 c6                	mov    %eax,%esi
   137ae:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   137b2:	48 8b 00             	mov    (%rax),%rax
   137b5:	49 89 d0             	mov    %rdx,%r8
   137b8:	b9 00 00 00 00       	mov    $0x0,%ecx
   137bd:	48 8d 15 5c a0 00 00 	lea    0xa05c(%rip),%rdx        # 1d820 <ompi_mpi_byte@@Base>
   137c4:	48 89 c7             	mov    %rax,%rdi
   137c7:	e8 84 66 ff ff       	call   9e50 <MPI_Bcast@plt>
   137cc:	c7 45 a4 00 00 00 00 	movl   $0x0,-0x5c(%rbp)
   137d3:	e9 c0 00 00 00       	jmp    13898 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x86c>
   137d8:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   137db:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   137df:	48 98                	cltq   
   137e1:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   137e8:	00 
   137e9:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   137ed:	48 8b 10             	mov    (%rax),%rdx
   137f0:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   137f3:	48 98                	cltq   
   137f5:	48 c1 e0 06          	shl    $0x6,%rax
   137f9:	48 01 d0             	add    %rdx,%rax
   137fc:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
   13800:	48 01 ca             	add    %rcx,%rdx
   13803:	48 89 50 08          	mov    %rdx,0x8(%rax)
   13807:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   1380a:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   1380e:	48 98                	cltq   
   13810:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13817:	00 
   13818:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   1381c:	48 8b 10             	mov    (%rax),%rdx
   1381f:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   13822:	48 98                	cltq   
   13824:	48 c1 e0 06          	shl    $0x6,%rax
   13828:	48 01 d0             	add    %rdx,%rax
   1382b:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   1382f:	48 01 ca             	add    %rcx,%rdx
   13832:	48 89 50 10          	mov    %rdx,0x10(%rax)
   13836:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   13839:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   1383d:	48 98                	cltq   
   1383f:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13846:	00 
   13847:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   1384b:	48 8b 10             	mov    (%rax),%rdx
   1384e:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   13851:	48 98                	cltq   
   13853:	48 c1 e0 06          	shl    $0x6,%rax
   13857:	48 01 d0             	add    %rdx,%rax
   1385a:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   1385e:	48 01 ca             	add    %rcx,%rdx
   13861:	48 89 50 18          	mov    %rdx,0x18(%rax)
   13865:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   13868:	0f af 45 b4          	imul   -0x4c(%rbp),%eax
   1386c:	48 98                	cltq   
   1386e:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13875:	00 
   13876:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   1387a:	48 8b 10             	mov    (%rax),%rdx
   1387d:	8b 45 a4             	mov    -0x5c(%rbp),%eax
   13880:	48 98                	cltq   
   13882:	48 c1 e0 06          	shl    $0x6,%rax
   13886:	48 01 d0             	add    %rdx,%rax
   13889:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   1388d:	48 01 ca             	add    %rcx,%rdx
   13890:	48 89 50 20          	mov    %rdx,0x20(%rax)
   13894:	83 45 a4 01          	addl   $0x1,-0x5c(%rbp)
   13898:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   1389c:	8b 00                	mov    (%rax),%eax
   1389e:	39 45 a4             	cmp    %eax,-0x5c(%rbp)
   138a1:	0f 8c 31 ff ff ff    	jl     137d8 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0x7ac>
   138a7:	48 8b 45 10          	mov    0x10(%rbp),%rax
   138ab:	48 8b 50 10          	mov    0x10(%rax),%rdx
   138af:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   138b3:	8b 00                	mov    (%rax),%eax
   138b5:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   138b9:	89 c6                	mov    %eax,%esi
   138bb:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   138bf:	48 8b 00             	mov    (%rax),%rax
   138c2:	48 8b 40 08          	mov    0x8(%rax),%rax
   138c6:	49 89 d0             	mov    %rdx,%r8
   138c9:	b9 00 00 00 00       	mov    $0x0,%ecx
   138ce:	48 8d 15 4b a3 00 00 	lea    0xa34b(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   138d5:	48 89 c7             	mov    %rax,%rdi
   138d8:	e8 73 65 ff ff       	call   9e50 <MPI_Bcast@plt>
   138dd:	48 8b 45 10          	mov    0x10(%rbp),%rax
   138e1:	48 8b 50 10          	mov    0x10(%rax),%rdx
   138e5:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   138e9:	8b 00                	mov    (%rax),%eax
   138eb:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   138ef:	89 c6                	mov    %eax,%esi
   138f1:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   138f5:	48 8b 00             	mov    (%rax),%rax
   138f8:	48 8b 40 10          	mov    0x10(%rax),%rax
   138fc:	49 89 d0             	mov    %rdx,%r8
   138ff:	b9 00 00 00 00       	mov    $0x0,%ecx
   13904:	48 8d 15 15 a3 00 00 	lea    0xa315(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   1390b:	48 89 c7             	mov    %rax,%rdi
   1390e:	e8 3d 65 ff ff       	call   9e50 <MPI_Bcast@plt>
   13913:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13917:	48 8b 50 10          	mov    0x10(%rax),%rdx
   1391b:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   1391f:	8b 00                	mov    (%rax),%eax
   13921:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13925:	89 c6                	mov    %eax,%esi
   13927:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   1392b:	48 8b 00             	mov    (%rax),%rax
   1392e:	48 8b 40 18          	mov    0x18(%rax),%rax
   13932:	49 89 d0             	mov    %rdx,%r8
   13935:	b9 00 00 00 00       	mov    $0x0,%ecx
   1393a:	48 8d 15 df a2 00 00 	lea    0xa2df(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13941:	48 89 c7             	mov    %rax,%rdi
   13944:	e8 07 65 ff ff       	call   9e50 <MPI_Bcast@plt>
   13949:	48 8b 45 10          	mov    0x10(%rbp),%rax
   1394d:	48 8b 50 10          	mov    0x10(%rax),%rdx
   13951:	48 8b 45 80          	mov    -0x80(%rbp),%rax
   13955:	8b 00                	mov    (%rax),%eax
   13957:	0f af 45 b4          	imul   -0x4c(%rbp),%eax
   1395b:	89 c6                	mov    %eax,%esi
   1395d:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   13961:	48 8b 00             	mov    (%rax),%rax
   13964:	48 8b 40 20          	mov    0x20(%rax),%rax
   13968:	49 89 d0             	mov    %rdx,%r8
   1396b:	b9 00 00 00 00       	mov    $0x0,%ecx
   13970:	48 8d 15 a9 a2 00 00 	lea    0xa2a9(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13977:	48 89 c7             	mov    %rax,%rdi
   1397a:	e8 d1 64 ff ff       	call   9e50 <MPI_Bcast@plt>
   1397f:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   13983:	8b 00                	mov    (%rax),%eax
   13985:	85 c0                	test   %eax,%eax
   13987:	0f 85 4c 03 00 00    	jne    13cd9 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0xcad>
   1398d:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13994:	48 8b 00             	mov    (%rax),%rax
   13997:	48 8b 40 08          	mov    0x8(%rax),%rax
   1399b:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   1399f:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   139a6:	48 8b 00             	mov    (%rax),%rax
   139a9:	48 8b 40 10          	mov    0x10(%rax),%rax
   139ad:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   139b1:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   139b8:	48 8b 00             	mov    (%rax),%rax
   139bb:	48 8b 40 18          	mov    0x18(%rax),%rax
   139bf:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   139c3:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   139ca:	48 8b 00             	mov    (%rax),%rax
   139cd:	48 8b 40 20          	mov    0x20(%rax),%rax
   139d1:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   139d5:	48 8b 45 10          	mov    0x10(%rbp),%rax
   139d9:	48 8b 70 28          	mov    0x28(%rax),%rsi
   139dd:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   139e4:	8b 00                	mov    (%rax),%eax
   139e6:	48 98                	cltq   
   139e8:	c1 e0 06             	shl    $0x6,%eax
   139eb:	41 89 c0             	mov    %eax,%r8d
   139ee:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   139f5:	48 8b 08             	mov    (%rax),%rcx
   139f8:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   139fe:	48 98                	cltq   
   13a00:	48 c1 e0 06          	shl    $0x6,%rax
   13a04:	48 89 c7             	mov    %rax,%rdi
   13a07:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   13a0b:	8b 00                	mov    (%rax),%eax
   13a0d:	4c 63 d8             	movslq %eax,%r11
   13a10:	48 89 f8             	mov    %rdi,%rax
   13a13:	ba 00 00 00 00       	mov    $0x0,%edx
   13a18:	49 f7 f3             	div    %r11
   13a1b:	89 c7                	mov    %eax,%edi
   13a1d:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   13a24:	56                   	push   %rsi
   13a25:	6a 00                	push   $0x0
   13a27:	4c 8d 0d f2 9d 00 00 	lea    0x9df2(%rip),%r9        # 1d820 <ompi_mpi_byte@@Base>
   13a2e:	48 8d 15 eb 9d 00 00 	lea    0x9deb(%rip),%rdx        # 1d820 <ompi_mpi_byte@@Base>
   13a35:	89 fe                	mov    %edi,%esi
   13a37:	48 89 c7             	mov    %rax,%rdi
   13a3a:	e8 11 67 ff ff       	call   a150 <MPI_Scatter@plt>
   13a3f:	48 83 c4 10          	add    $0x10,%rsp
   13a43:	c7 45 a8 00 00 00 00 	movl   $0x0,-0x58(%rbp)
   13a4a:	e9 cc 00 00 00       	jmp    13b1b <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0xaef>
   13a4f:	8b 45 a8             	mov    -0x58(%rbp),%eax
   13a52:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13a56:	48 98                	cltq   
   13a58:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13a5f:	00 
   13a60:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13a67:	48 8b 10             	mov    (%rax),%rdx
   13a6a:	8b 45 a8             	mov    -0x58(%rbp),%eax
   13a6d:	48 98                	cltq   
   13a6f:	48 c1 e0 06          	shl    $0x6,%rax
   13a73:	48 01 d0             	add    %rdx,%rax
   13a76:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
   13a7a:	48 01 ca             	add    %rcx,%rdx
   13a7d:	48 89 50 08          	mov    %rdx,0x8(%rax)
   13a81:	8b 45 a8             	mov    -0x58(%rbp),%eax
   13a84:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13a88:	48 98                	cltq   
   13a8a:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13a91:	00 
   13a92:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13a99:	48 8b 10             	mov    (%rax),%rdx
   13a9c:	8b 45 a8             	mov    -0x58(%rbp),%eax
   13a9f:	48 98                	cltq   
   13aa1:	48 c1 e0 06          	shl    $0x6,%rax
   13aa5:	48 01 d0             	add    %rdx,%rax
   13aa8:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   13aac:	48 01 ca             	add    %rcx,%rdx
   13aaf:	48 89 50 10          	mov    %rdx,0x10(%rax)
   13ab3:	8b 45 a8             	mov    -0x58(%rbp),%eax
   13ab6:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13aba:	48 98                	cltq   
   13abc:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13ac3:	00 
   13ac4:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13acb:	48 8b 10             	mov    (%rax),%rdx
   13ace:	8b 45 a8             	mov    -0x58(%rbp),%eax
   13ad1:	48 98                	cltq   
   13ad3:	48 c1 e0 06          	shl    $0x6,%rax
   13ad7:	48 01 d0             	add    %rdx,%rax
   13ada:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   13ade:	48 01 ca             	add    %rcx,%rdx
   13ae1:	48 89 50 18          	mov    %rdx,0x18(%rax)
   13ae5:	8b 45 a8             	mov    -0x58(%rbp),%eax
   13ae8:	0f af 45 b4          	imul   -0x4c(%rbp),%eax
   13aec:	48 98                	cltq   
   13aee:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13af5:	00 
   13af6:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13afd:	48 8b 10             	mov    (%rax),%rdx
   13b00:	8b 45 a8             	mov    -0x58(%rbp),%eax
   13b03:	48 98                	cltq   
   13b05:	48 c1 e0 06          	shl    $0x6,%rax
   13b09:	48 01 d0             	add    %rdx,%rax
   13b0c:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   13b10:	48 01 ca             	add    %rcx,%rdx
   13b13:	48 89 50 20          	mov    %rdx,0x20(%rax)
   13b17:	83 45 a8 01          	addl   $0x1,-0x58(%rbp)
   13b1b:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   13b22:	8b 00                	mov    (%rax),%eax
   13b24:	39 45 a8             	cmp    %eax,-0x58(%rbp)
   13b27:	0f 8c 22 ff ff ff    	jl     13a4f <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0xa23>
   13b2d:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13b31:	48 8b 70 28          	mov    0x28(%rax),%rsi
   13b35:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   13b3c:	8b 00                	mov    (%rax),%eax
   13b3e:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13b42:	89 c7                	mov    %eax,%edi
   13b44:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13b4b:	48 8b 00             	mov    (%rax),%rax
   13b4e:	48 8b 48 08          	mov    0x8(%rax),%rcx
   13b52:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   13b58:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13b5c:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
   13b60:	44 8b 1a             	mov    (%rdx),%r11d
   13b63:	99                   	cltd   
   13b64:	41 f7 fb             	idiv   %r11d
   13b67:	41 89 c2             	mov    %eax,%r10d
   13b6a:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   13b71:	48 8b 40 08          	mov    0x8(%rax),%rax
   13b75:	56                   	push   %rsi
   13b76:	6a 00                	push   $0x0
   13b78:	4c 8d 0d a1 a0 00 00 	lea    0xa0a1(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   13b7f:	41 89 f8             	mov    %edi,%r8d
   13b82:	48 8d 15 97 a0 00 00 	lea    0xa097(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13b89:	44 89 d6             	mov    %r10d,%esi
   13b8c:	48 89 c7             	mov    %rax,%rdi
   13b8f:	e8 bc 65 ff ff       	call   a150 <MPI_Scatter@plt>
   13b94:	48 83 c4 10          	add    $0x10,%rsp
   13b98:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13b9c:	48 8b 70 28          	mov    0x28(%rax),%rsi
   13ba0:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   13ba7:	8b 00                	mov    (%rax),%eax
   13ba9:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13bad:	89 c7                	mov    %eax,%edi
   13baf:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13bb6:	48 8b 00             	mov    (%rax),%rax
   13bb9:	48 8b 48 10          	mov    0x10(%rax),%rcx
   13bbd:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   13bc3:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13bc7:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
   13bcb:	44 8b 1a             	mov    (%rdx),%r11d
   13bce:	99                   	cltd   
   13bcf:	41 f7 fb             	idiv   %r11d
   13bd2:	41 89 c2             	mov    %eax,%r10d
   13bd5:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   13bdc:	48 8b 40 10          	mov    0x10(%rax),%rax
   13be0:	56                   	push   %rsi
   13be1:	6a 00                	push   $0x0
   13be3:	4c 8d 0d 36 a0 00 00 	lea    0xa036(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   13bea:	41 89 f8             	mov    %edi,%r8d
   13bed:	48 8d 15 2c a0 00 00 	lea    0xa02c(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13bf4:	44 89 d6             	mov    %r10d,%esi
   13bf7:	48 89 c7             	mov    %rax,%rdi
   13bfa:	e8 51 65 ff ff       	call   a150 <MPI_Scatter@plt>
   13bff:	48 83 c4 10          	add    $0x10,%rsp
   13c03:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13c07:	48 8b 70 28          	mov    0x28(%rax),%rsi
   13c0b:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   13c12:	8b 00                	mov    (%rax),%eax
   13c14:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13c18:	89 c7                	mov    %eax,%edi
   13c1a:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13c21:	48 8b 00             	mov    (%rax),%rax
   13c24:	48 8b 48 18          	mov    0x18(%rax),%rcx
   13c28:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   13c2e:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13c32:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
   13c36:	44 8b 12             	mov    (%rdx),%r10d
   13c39:	99                   	cltd   
   13c3a:	41 f7 fa             	idiv   %r10d
   13c3d:	41 89 c2             	mov    %eax,%r10d
   13c40:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   13c47:	48 8b 40 18          	mov    0x18(%rax),%rax
   13c4b:	56                   	push   %rsi
   13c4c:	6a 00                	push   $0x0
   13c4e:	4c 8d 0d cb 9f 00 00 	lea    0x9fcb(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   13c55:	41 89 f8             	mov    %edi,%r8d
   13c58:	48 8d 15 c1 9f 00 00 	lea    0x9fc1(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13c5f:	44 89 d6             	mov    %r10d,%esi
   13c62:	48 89 c7             	mov    %rax,%rdi
   13c65:	e8 e6 64 ff ff       	call   a150 <MPI_Scatter@plt>
   13c6a:	48 83 c4 10          	add    $0x10,%rsp
   13c6e:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13c72:	48 8b 70 28          	mov    0x28(%rax),%rsi
   13c76:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   13c7d:	8b 00                	mov    (%rax),%eax
   13c7f:	0f af 45 b4          	imul   -0x4c(%rbp),%eax
   13c83:	89 c7                	mov    %eax,%edi
   13c85:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13c8c:	48 8b 00             	mov    (%rax),%rax
   13c8f:	48 8b 48 20          	mov    0x20(%rax),%rcx
   13c93:	8b 85 64 ff ff ff    	mov    -0x9c(%rbp),%eax
   13c99:	0f af 45 b4          	imul   -0x4c(%rbp),%eax
   13c9d:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
   13ca1:	44 8b 1a             	mov    (%rdx),%r11d
   13ca4:	99                   	cltd   
   13ca5:	41 f7 fb             	idiv   %r11d
   13ca8:	41 89 c2             	mov    %eax,%r10d
   13cab:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   13cb2:	48 8b 40 20          	mov    0x20(%rax),%rax
   13cb6:	56                   	push   %rsi
   13cb7:	6a 00                	push   $0x0
   13cb9:	4c 8d 0d 60 9f 00 00 	lea    0x9f60(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   13cc0:	41 89 f8             	mov    %edi,%r8d
   13cc3:	48 8d 15 56 9f 00 00 	lea    0x9f56(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13cca:	44 89 d6             	mov    %r10d,%esi
   13ccd:	48 89 c7             	mov    %rax,%rdi
   13cd0:	e8 7b 64 ff ff       	call   a150 <MPI_Scatter@plt>
   13cd5:	48 83 c4 10          	add    $0x10,%rsp
   13cd9:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13ce0:	48 8b 00             	mov    (%rax),%rax
   13ce3:	48 8b 40 08          	mov    0x8(%rax),%rax
   13ce7:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   13ceb:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13cf2:	48 8b 00             	mov    (%rax),%rax
   13cf5:	48 8b 40 10          	mov    0x10(%rax),%rax
   13cf9:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   13cfd:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13d04:	48 8b 00             	mov    (%rax),%rax
   13d07:	48 8b 40 18          	mov    0x18(%rax),%rax
   13d0b:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   13d0f:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13d16:	48 8b 00             	mov    (%rax),%rax
   13d19:	48 8b 40 20          	mov    0x20(%rax),%rax
   13d1d:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   13d21:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13d25:	48 8b 50 18          	mov    0x18(%rax),%rdx
   13d29:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   13d30:	8b 00                	mov    (%rax),%eax
   13d32:	48 98                	cltq   
   13d34:	c1 e0 06             	shl    $0x6,%eax
   13d37:	89 c6                	mov    %eax,%esi
   13d39:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13d40:	48 8b 00             	mov    (%rax),%rax
   13d43:	49 89 d0             	mov    %rdx,%r8
   13d46:	b9 00 00 00 00       	mov    $0x0,%ecx
   13d4b:	48 8d 15 ce 9a 00 00 	lea    0x9ace(%rip),%rdx        # 1d820 <ompi_mpi_byte@@Base>
   13d52:	48 89 c7             	mov    %rax,%rdi
   13d55:	e8 f6 60 ff ff       	call   9e50 <MPI_Bcast@plt>
   13d5a:	c7 45 ac 00 00 00 00 	movl   $0x0,-0x54(%rbp)
   13d61:	e9 cc 00 00 00       	jmp    13e32 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0xe06>
   13d66:	8b 45 ac             	mov    -0x54(%rbp),%eax
   13d69:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13d6d:	48 98                	cltq   
   13d6f:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13d76:	00 
   13d77:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13d7e:	48 8b 10             	mov    (%rax),%rdx
   13d81:	8b 45 ac             	mov    -0x54(%rbp),%eax
   13d84:	48 98                	cltq   
   13d86:	48 c1 e0 06          	shl    $0x6,%rax
   13d8a:	48 01 d0             	add    %rdx,%rax
   13d8d:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
   13d91:	48 01 ca             	add    %rcx,%rdx
   13d94:	48 89 50 08          	mov    %rdx,0x8(%rax)
   13d98:	8b 45 ac             	mov    -0x54(%rbp),%eax
   13d9b:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13d9f:	48 98                	cltq   
   13da1:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13da8:	00 
   13da9:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13db0:	48 8b 10             	mov    (%rax),%rdx
   13db3:	8b 45 ac             	mov    -0x54(%rbp),%eax
   13db6:	48 98                	cltq   
   13db8:	48 c1 e0 06          	shl    $0x6,%rax
   13dbc:	48 01 d0             	add    %rdx,%rax
   13dbf:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   13dc3:	48 01 ca             	add    %rcx,%rdx
   13dc6:	48 89 50 10          	mov    %rdx,0x10(%rax)
   13dca:	8b 45 ac             	mov    -0x54(%rbp),%eax
   13dcd:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13dd1:	48 98                	cltq   
   13dd3:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13dda:	00 
   13ddb:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13de2:	48 8b 10             	mov    (%rax),%rdx
   13de5:	8b 45 ac             	mov    -0x54(%rbp),%eax
   13de8:	48 98                	cltq   
   13dea:	48 c1 e0 06          	shl    $0x6,%rax
   13dee:	48 01 d0             	add    %rdx,%rax
   13df1:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   13df5:	48 01 ca             	add    %rcx,%rdx
   13df8:	48 89 50 18          	mov    %rdx,0x18(%rax)
   13dfc:	8b 45 ac             	mov    -0x54(%rbp),%eax
   13dff:	0f af 45 b4          	imul   -0x4c(%rbp),%eax
   13e03:	48 98                	cltq   
   13e05:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   13e0c:	00 
   13e0d:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13e14:	48 8b 10             	mov    (%rax),%rdx
   13e17:	8b 45 ac             	mov    -0x54(%rbp),%eax
   13e1a:	48 98                	cltq   
   13e1c:	48 c1 e0 06          	shl    $0x6,%rax
   13e20:	48 01 d0             	add    %rdx,%rax
   13e23:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   13e27:	48 01 ca             	add    %rcx,%rdx
   13e2a:	48 89 50 20          	mov    %rdx,0x20(%rax)
   13e2e:	83 45 ac 01          	addl   $0x1,-0x54(%rbp)
   13e32:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   13e39:	8b 00                	mov    (%rax),%eax
   13e3b:	39 45 ac             	cmp    %eax,-0x54(%rbp)
   13e3e:	0f 8c 22 ff ff ff    	jl     13d66 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0xd3a>
   13e44:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13e48:	48 8b 50 18          	mov    0x18(%rax),%rdx
   13e4c:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   13e53:	8b 00                	mov    (%rax),%eax
   13e55:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13e59:	89 c6                	mov    %eax,%esi
   13e5b:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13e62:	48 8b 00             	mov    (%rax),%rax
   13e65:	48 8b 40 08          	mov    0x8(%rax),%rax
   13e69:	49 89 d0             	mov    %rdx,%r8
   13e6c:	b9 00 00 00 00       	mov    $0x0,%ecx
   13e71:	48 8d 15 a8 9d 00 00 	lea    0x9da8(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13e78:	48 89 c7             	mov    %rax,%rdi
   13e7b:	e8 d0 5f ff ff       	call   9e50 <MPI_Bcast@plt>
   13e80:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13e84:	48 8b 50 18          	mov    0x18(%rax),%rdx
   13e88:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   13e8f:	8b 00                	mov    (%rax),%eax
   13e91:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13e95:	89 c6                	mov    %eax,%esi
   13e97:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13e9e:	48 8b 00             	mov    (%rax),%rax
   13ea1:	48 8b 40 10          	mov    0x10(%rax),%rax
   13ea5:	49 89 d0             	mov    %rdx,%r8
   13ea8:	b9 00 00 00 00       	mov    $0x0,%ecx
   13ead:	48 8d 15 6c 9d 00 00 	lea    0x9d6c(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13eb4:	48 89 c7             	mov    %rax,%rdi
   13eb7:	e8 94 5f ff ff       	call   9e50 <MPI_Bcast@plt>
   13ebc:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13ec0:	48 8b 50 18          	mov    0x18(%rax),%rdx
   13ec4:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   13ecb:	8b 00                	mov    (%rax),%eax
   13ecd:	0f af 45 b0          	imul   -0x50(%rbp),%eax
   13ed1:	89 c6                	mov    %eax,%esi
   13ed3:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13eda:	48 8b 00             	mov    (%rax),%rax
   13edd:	48 8b 40 18          	mov    0x18(%rax),%rax
   13ee1:	49 89 d0             	mov    %rdx,%r8
   13ee4:	b9 00 00 00 00       	mov    $0x0,%ecx
   13ee9:	48 8d 15 30 9d 00 00 	lea    0x9d30(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13ef0:	48 89 c7             	mov    %rax,%rdi
   13ef3:	e8 58 5f ff ff       	call   9e50 <MPI_Bcast@plt>
   13ef8:	48 8b 45 10          	mov    0x10(%rbp),%rax
   13efc:	48 8b 50 18          	mov    0x18(%rax),%rdx
   13f00:	48 8b 85 70 ff ff ff 	mov    -0x90(%rbp),%rax
   13f07:	8b 00                	mov    (%rax),%eax
   13f09:	0f af 45 b4          	imul   -0x4c(%rbp),%eax
   13f0d:	89 c6                	mov    %eax,%esi
   13f0f:	48 8b 85 78 ff ff ff 	mov    -0x88(%rbp),%rax
   13f16:	48 8b 00             	mov    (%rax),%rax
   13f19:	48 8b 40 20          	mov    0x20(%rax),%rax
   13f1d:	49 89 d0             	mov    %rdx,%r8
   13f20:	b9 00 00 00 00       	mov    $0x0,%ecx
   13f25:	48 8d 15 f4 9c 00 00 	lea    0x9cf4(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   13f2c:	48 89 c7             	mov    %rax,%rdi
   13f2f:	e8 1c 5f ff ff       	call   9e50 <MPI_Bcast@plt>
   13f34:	80 7d 9f 00          	cmpb   $0x0,-0x61(%rbp)
   13f38:	74 1a                	je     13f54 <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t+0xf28>
   13f3a:	48 8b 85 68 ff ff ff 	mov    -0x98(%rbp),%rax
   13f41:	48 89 c7             	mov    %rax,%rdi
   13f44:	e8 7f dc ff ff       	call   11bc8 <_Z14free_particlesP10particle_t>
   13f49:	48 c7 85 68 ff ff ff 	movq   $0x0,-0x98(%rbp)
   13f50:	00 00 00 00 
   13f54:	90                   	nop
   13f55:	c9                   	leave  
   13f56:	c3                   	ret    

0000000000013f57 <_Z7scatterPP17chunk_particles_tS1_S0_P10topology_t>:
   13f57:	f3 0f 1e fa          	endbr64 
   13f5b:	55                   	push   %rbp
   13f5c:	48 89 e5             	mov    %rsp,%rbp
   13f5f:	48 83 ec 50          	sub    $0x50,%rsp
   13f63:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
   13f67:	48 89 75 c0          	mov    %rsi,-0x40(%rbp)
   13f6b:	48 89 55 b8          	mov    %rdx,-0x48(%rbp)
   13f6f:	48 89 4d b0          	mov    %rcx,-0x50(%rbp)
   13f73:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   13f7a:	00 00 
   13f7c:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   13f80:	31 c0                	xor    %eax,%eax
   13f82:	bf 10 00 00 00       	mov    $0x10,%edi
   13f87:	e8 94 66 ff ff       	call   a620 <malloc@plt>
   13f8c:	48 89 c2             	mov    %rax,%rdx
   13f8f:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   13f93:	48 89 10             	mov    %rdx,(%rax)
   13f96:	bf 10 00 00 00       	mov    $0x10,%edi
   13f9b:	e8 80 66 ff ff       	call   a620 <malloc@plt>
   13fa0:	48 89 c2             	mov    %rax,%rdx
   13fa3:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   13fa7:	48 89 10             	mov    %rdx,(%rax)
   13faa:	48 c7 45 e0 00 00 00 	movq   $0x0,-0x20(%rbp)
   13fb1:	00 
   13fb2:	48 c7 45 e8 00 00 00 	movq   $0x0,-0x18(%rbp)
   13fb9:	00 
   13fba:	48 c7 45 f0 00 00 00 	movq   $0x0,-0x10(%rbp)
   13fc1:	00 
   13fc2:	c7 45 d4 ff ff ff ff 	movl   $0xffffffff,-0x2c(%rbp)
   13fc9:	c7 45 d8 ff ff ff ff 	movl   $0xffffffff,-0x28(%rbp)
   13fd0:	48 83 7d b8 00       	cmpq   $0x0,-0x48(%rbp)
   13fd5:	75 11                	jne    13fe8 <_Z7scatterPP17chunk_particles_tS1_S0_P10topology_t+0x91>
   13fd7:	48 c7 45 f0 00 00 00 	movq   $0x0,-0x10(%rbp)
   13fde:	00 
   13fdf:	c7 45 dc ff ff ff ff 	movl   $0xffffffff,-0x24(%rbp)
   13fe6:	eb 15                	jmp    13ffd <_Z7scatterPP17chunk_particles_tS1_S0_P10topology_t+0xa6>
   13fe8:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   13fec:	48 8b 00             	mov    (%rax),%rax
   13fef:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   13ff3:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   13ff7:	8b 40 08             	mov    0x8(%rax),%eax
   13ffa:	89 45 dc             	mov    %eax,-0x24(%rbp)
   13ffd:	44 8b 45 dc          	mov    -0x24(%rbp),%r8d
   14001:	48 8b 7d f0          	mov    -0x10(%rbp),%rdi
   14005:	48 8d 4d d8          	lea    -0x28(%rbp),%rcx
   14009:	48 8d 55 e8          	lea    -0x18(%rbp),%rdx
   1400d:	48 8d 75 d4          	lea    -0x2c(%rbp),%rsi
   14011:	48 8d 45 e0          	lea    -0x20(%rbp),%rax
   14015:	48 83 ec 08          	sub    $0x8,%rsp
   14019:	ff 75 b0             	push   -0x50(%rbp)
   1401c:	45 89 c1             	mov    %r8d,%r9d
   1401f:	49 89 f8             	mov    %rdi,%r8
   14022:	48 89 c7             	mov    %rax,%rdi
   14025:	e8 02 f0 ff ff       	call   1302c <_Z7scatterPP10particle_tPiS1_S2_S0_iP10topology_t>
   1402a:	48 83 c4 10          	add    $0x10,%rsp
   1402e:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14032:	48 8b 00             	mov    (%rax),%rax
   14035:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
   14039:	48 89 10             	mov    %rdx,(%rax)
   1403c:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14040:	48 8b 00             	mov    (%rax),%rax
   14043:	8b 55 d4             	mov    -0x2c(%rbp),%edx
   14046:	89 50 08             	mov    %edx,0x8(%rax)
   14049:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   1404d:	48 8b 00             	mov    (%rax),%rax
   14050:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   14054:	48 89 10             	mov    %rdx,(%rax)
   14057:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   1405b:	48 8b 00             	mov    (%rax),%rax
   1405e:	8b 55 d8             	mov    -0x28(%rbp),%edx
   14061:	89 50 08             	mov    %edx,0x8(%rax)
   14064:	90                   	nop
   14065:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14069:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   14070:	00 00 
   14072:	74 05                	je     14079 <_Z7scatterPP17chunk_particles_tS1_S0_P10topology_t+0x122>
   14074:	e8 57 64 ff ff       	call   a4d0 <__stack_chk_fail@plt>
   14079:	c9                   	leave  
   1407a:	c3                   	ret    

000000000001407b <_Z6reduceP10particle_tiS0_iP10topology_t>:
   1407b:	f3 0f 1e fa          	endbr64 
   1407f:	55                   	push   %rbp
   14080:	48 89 e5             	mov    %rsp,%rbp
   14083:	48 83 ec 40          	sub    $0x40,%rsp
   14087:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   1408b:	89 75 d4             	mov    %esi,-0x2c(%rbp)
   1408e:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
   14092:	89 4d d0             	mov    %ecx,-0x30(%rbp)
   14095:	4c 89 45 c0          	mov    %r8,-0x40(%rbp)
   14099:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   1409d:	8b 40 30             	mov    0x30(%rax),%eax
   140a0:	0f af 45 d4          	imul   -0x2c(%rbp),%eax
   140a4:	48 98                	cltq   
   140a6:	48 c1 e0 03          	shl    $0x3,%rax
   140aa:	48 89 c7             	mov    %rax,%rdi
   140ad:	e8 6e 65 ff ff       	call   a620 <malloc@plt>
   140b2:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   140b6:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   140ba:	48 8b 48 18          	mov    0x18(%rax),%rcx
   140be:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   140c2:	8b 40 30             	mov    0x30(%rax),%eax
   140c5:	0f af 45 d4          	imul   -0x2c(%rbp),%eax
   140c9:	89 c2                	mov    %eax,%edx
   140cb:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   140cf:	48 8b 40 18          	mov    0x18(%rax),%rax
   140d3:	48 8b 75 f0          	mov    -0x10(%rbp),%rsi
   140d7:	49 89 c9             	mov    %rcx,%r9
   140da:	4c 8d 05 3f 8f 00 00 	lea    0x8f3f(%rip),%r8        # 1d020 <ompi_mpi_op_sum@@Base>
   140e1:	48 8d 0d 38 9b 00 00 	lea    0x9b38(%rip),%rcx        # 1dc20 <ompi_mpi_double@@Base>
   140e8:	48 89 c7             	mov    %rax,%rdi
   140eb:	e8 60 64 ff ff       	call   a550 <MPI_Allreduce@plt>
   140f0:	c7 45 e0 00 00 00 00 	movl   $0x0,-0x20(%rbp)
   140f7:	eb 6f                	jmp    14168 <_Z6reduceP10particle_tiS0_iP10topology_t+0xed>
   140f9:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
   14100:	eb 56                	jmp    14158 <_Z6reduceP10particle_tiS0_iP10topology_t+0xdd>
   14102:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   14106:	8b 40 30             	mov    0x30(%rax),%eax
   14109:	0f af 45 e0          	imul   -0x20(%rbp),%eax
   1410d:	48 63 d0             	movslq %eax,%rdx
   14110:	8b 45 e4             	mov    -0x1c(%rbp),%eax
   14113:	48 98                	cltq   
   14115:	48 01 d0             	add    %rdx,%rax
   14118:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   1411f:	00 
   14120:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   14124:	48 01 c2             	add    %rax,%rdx
   14127:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   1412b:	48 8b 48 18          	mov    0x18(%rax),%rcx
   1412f:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   14133:	8b 40 30             	mov    0x30(%rax),%eax
   14136:	0f af 45 e0          	imul   -0x20(%rbp),%eax
   1413a:	48 63 f0             	movslq %eax,%rsi
   1413d:	8b 45 e4             	mov    -0x1c(%rbp),%eax
   14140:	48 98                	cltq   
   14142:	48 01 f0             	add    %rsi,%rax
   14145:	48 c1 e0 03          	shl    $0x3,%rax
   14149:	48 01 c8             	add    %rcx,%rax
   1414c:	f2 0f 10 02          	movsd  (%rdx),%xmm0
   14150:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   14154:	83 45 e4 01          	addl   $0x1,-0x1c(%rbp)
   14158:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   1415c:	8b 40 30             	mov    0x30(%rax),%eax
   1415f:	39 45 e4             	cmp    %eax,-0x1c(%rbp)
   14162:	7c 9e                	jl     14102 <_Z6reduceP10particle_tiS0_iP10topology_t+0x87>
   14164:	83 45 e0 01          	addl   $0x1,-0x20(%rbp)
   14168:	8b 45 e0             	mov    -0x20(%rbp),%eax
   1416b:	3b 45 d4             	cmp    -0x2c(%rbp),%eax
   1416e:	7c 89                	jl     140f9 <_Z6reduceP10particle_tiS0_iP10topology_t+0x7e>
   14170:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14174:	8b 40 30             	mov    0x30(%rax),%eax
   14177:	0f af 45 d0          	imul   -0x30(%rbp),%eax
   1417b:	48 98                	cltq   
   1417d:	48 c1 e0 03          	shl    $0x3,%rax
   14181:	48 89 c7             	mov    %rax,%rdi
   14184:	e8 97 64 ff ff       	call   a620 <malloc@plt>
   14189:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   1418d:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   14191:	48 8b 48 10          	mov    0x10(%rax),%rcx
   14195:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14199:	8b 40 30             	mov    0x30(%rax),%eax
   1419c:	0f af 45 d0          	imul   -0x30(%rbp),%eax
   141a0:	89 c2                	mov    %eax,%edx
   141a2:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   141a6:	48 8b 40 18          	mov    0x18(%rax),%rax
   141aa:	48 8b 75 f8          	mov    -0x8(%rbp),%rsi
   141ae:	49 89 c9             	mov    %rcx,%r9
   141b1:	4c 8d 05 68 8e 00 00 	lea    0x8e68(%rip),%r8        # 1d020 <ompi_mpi_op_sum@@Base>
   141b8:	48 8d 0d 61 9a 00 00 	lea    0x9a61(%rip),%rcx        # 1dc20 <ompi_mpi_double@@Base>
   141bf:	48 89 c7             	mov    %rax,%rdi
   141c2:	e8 89 63 ff ff       	call   a550 <MPI_Allreduce@plt>
   141c7:	c7 45 e8 00 00 00 00 	movl   $0x0,-0x18(%rbp)
   141ce:	eb 6f                	jmp    1423f <_Z6reduceP10particle_tiS0_iP10topology_t+0x1c4>
   141d0:	c7 45 ec 00 00 00 00 	movl   $0x0,-0x14(%rbp)
   141d7:	eb 56                	jmp    1422f <_Z6reduceP10particle_tiS0_iP10topology_t+0x1b4>
   141d9:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   141dd:	8b 40 30             	mov    0x30(%rax),%eax
   141e0:	0f af 45 e8          	imul   -0x18(%rbp),%eax
   141e4:	48 63 d0             	movslq %eax,%rdx
   141e7:	8b 45 ec             	mov    -0x14(%rbp),%eax
   141ea:	48 98                	cltq   
   141ec:	48 01 d0             	add    %rdx,%rax
   141ef:	48 8d 14 c5 00 00 00 	lea    0x0(,%rax,8),%rdx
   141f6:	00 
   141f7:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   141fb:	48 01 c2             	add    %rax,%rdx
   141fe:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14202:	48 8b 48 18          	mov    0x18(%rax),%rcx
   14206:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   1420a:	8b 40 30             	mov    0x30(%rax),%eax
   1420d:	0f af 45 e8          	imul   -0x18(%rbp),%eax
   14211:	48 63 f0             	movslq %eax,%rsi
   14214:	8b 45 ec             	mov    -0x14(%rbp),%eax
   14217:	48 98                	cltq   
   14219:	48 01 f0             	add    %rsi,%rax
   1421c:	48 c1 e0 03          	shl    $0x3,%rax
   14220:	48 01 c8             	add    %rcx,%rax
   14223:	f2 0f 10 02          	movsd  (%rdx),%xmm0
   14227:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   1422b:	83 45 ec 01          	addl   $0x1,-0x14(%rbp)
   1422f:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14233:	8b 40 30             	mov    0x30(%rax),%eax
   14236:	39 45 ec             	cmp    %eax,-0x14(%rbp)
   14239:	7c 9e                	jl     141d9 <_Z6reduceP10particle_tiS0_iP10topology_t+0x15e>
   1423b:	83 45 e8 01          	addl   $0x1,-0x18(%rbp)
   1423f:	8b 45 e8             	mov    -0x18(%rbp),%eax
   14242:	3b 45 d0             	cmp    -0x30(%rbp),%eax
   14245:	7c 89                	jl     141d0 <_Z6reduceP10particle_tiS0_iP10topology_t+0x155>
   14247:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   1424b:	48 89 c7             	mov    %rax,%rdi
   1424e:	e8 9d 62 ff ff       	call   a4f0 <free@plt>
   14253:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14257:	48 89 c7             	mov    %rax,%rdi
   1425a:	e8 91 62 ff ff       	call   a4f0 <free@plt>
   1425f:	90                   	nop
   14260:	c9                   	leave  
   14261:	c3                   	ret    

0000000000014262 <_Z6reduceP17chunk_particles_tS0_P10topology_t>:
   14262:	f3 0f 1e fa          	endbr64 
   14266:	55                   	push   %rbp
   14267:	48 89 e5             	mov    %rsp,%rbp
   1426a:	48 83 ec 20          	sub    $0x20,%rsp
   1426e:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   14272:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   14276:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
   1427a:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   1427e:	8b 48 08             	mov    0x8(%rax),%ecx
   14281:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   14285:	48 8b 10             	mov    (%rax),%rdx
   14288:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1428c:	8b 70 08             	mov    0x8(%rax),%esi
   1428f:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14293:	48 8b 00             	mov    (%rax),%rax
   14296:	48 8b 7d e8          	mov    -0x18(%rbp),%rdi
   1429a:	49 89 f8             	mov    %rdi,%r8
   1429d:	48 89 c7             	mov    %rax,%rdi
   142a0:	e8 d6 fd ff ff       	call   1407b <_Z6reduceP10particle_tiS0_iP10topology_t>
   142a5:	90                   	nop
   142a6:	c9                   	leave  
   142a7:	c3                   	ret    

00000000000142a8 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t>:
   142a8:	f3 0f 1e fa          	endbr64 
   142ac:	55                   	push   %rbp
   142ad:	48 89 e5             	mov    %rsp,%rbp
   142b0:	48 83 c4 80          	add    $0xffffffffffffff80,%rsp
   142b4:	48 89 7d a8          	mov    %rdi,-0x58(%rbp)
   142b8:	89 75 a4             	mov    %esi,-0x5c(%rbp)
   142bb:	48 89 55 98          	mov    %rdx,-0x68(%rbp)
   142bf:	89 4d a0             	mov    %ecx,-0x60(%rbp)
   142c2:	4c 89 45 90          	mov    %r8,-0x70(%rbp)
   142c6:	4c 89 4d 88          	mov    %r9,-0x78(%rbp)
   142ca:	48 8b 45 10          	mov    0x10(%rbp),%rax
   142ce:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   142d2:	48 8b 45 10          	mov    0x10(%rbp),%rax
   142d6:	48 83 c0 04          	add    $0x4,%rax
   142da:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   142de:	48 8b 45 10          	mov    0x10(%rbp),%rax
   142e2:	48 83 c0 08          	add    $0x8,%rax
   142e6:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   142ea:	48 8b 45 10          	mov    0x10(%rbp),%rax
   142ee:	48 83 c0 0c          	add    $0xc,%rax
   142f2:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   142f6:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   142fa:	8b 00                	mov    (%rax),%eax
   142fc:	85 c0                	test   %eax,%eax
   142fe:	0f 85 82 00 00 00    	jne    14386 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0xde>
   14304:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   14308:	8b 00                	mov    (%rax),%eax
   1430a:	85 c0                	test   %eax,%eax
   1430c:	75 78                	jne    14386 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0xde>
   1430e:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14312:	8b 10                	mov    (%rax),%edx
   14314:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   14318:	8b 00                	mov    (%rax),%eax
   1431a:	0f af 45 a4          	imul   -0x5c(%rbp),%eax
   1431e:	39 c2                	cmp    %eax,%edx
   14320:	74 28                	je     1434a <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0xa2>
   14322:	48 8d 05 57 2b 00 00 	lea    0x2b57(%rip),%rax        # 16e80 <_ZN3MPIL11LOCK_SHAREDE+0x238>
   14329:	48 89 c1             	mov    %rax,%rcx
   1432c:	ba 37 02 00 00       	mov    $0x237,%edx
   14331:	48 8d 05 48 29 00 00 	lea    0x2948(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   14338:	48 89 c6             	mov    %rax,%rsi
   1433b:	48 8d 05 96 2b 00 00 	lea    0x2b96(%rip),%rax        # 16ed8 <_ZN3MPIL11LOCK_SHAREDE+0x290>
   14342:	48 89 c7             	mov    %rax,%rdi
   14345:	e8 e6 5d ff ff       	call   a130 <__assert_fail@plt>
   1434a:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   1434e:	8b 10                	mov    (%rax),%edx
   14350:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14354:	8b 00                	mov    (%rax),%eax
   14356:	0f af 45 a0          	imul   -0x60(%rbp),%eax
   1435a:	39 c2                	cmp    %eax,%edx
   1435c:	74 28                	je     14386 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0xde>
   1435e:	48 8d 05 1b 2b 00 00 	lea    0x2b1b(%rip),%rax        # 16e80 <_ZN3MPIL11LOCK_SHAREDE+0x238>
   14365:	48 89 c1             	mov    %rax,%rcx
   14368:	ba 38 02 00 00       	mov    $0x238,%edx
   1436d:	48 8d 05 0c 29 00 00 	lea    0x290c(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   14374:	48 89 c6             	mov    %rax,%rsi
   14377:	48 8d 05 8a 2b 00 00 	lea    0x2b8a(%rip),%rax        # 16f08 <_ZN3MPIL11LOCK_SHAREDE+0x2c0>
   1437e:	48 89 c7             	mov    %rax,%rdi
   14381:	e8 aa 5d ff ff       	call   a130 <__assert_fail@plt>
   14386:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   1438a:	8b 00                	mov    (%rax),%eax
   1438c:	85 c0                	test   %eax,%eax
   1438e:	0f 85 8e 03 00 00    	jne    14722 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x47a>
   14394:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   14398:	8b 00                	mov    (%rax),%eax
   1439a:	85 c0                	test   %eax,%eax
   1439c:	0f 85 2e 02 00 00    	jne    145d0 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x328>
   143a2:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   143a6:	48 8b 00             	mov    (%rax),%rax
   143a9:	48 8b 40 08          	mov    0x8(%rax),%rax
   143ad:	48 89 45 c0          	mov    %rax,-0x40(%rbp)
   143b1:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   143b5:	48 8b 00             	mov    (%rax),%rax
   143b8:	48 8b 40 10          	mov    0x10(%rax),%rax
   143bc:	48 89 45 c8          	mov    %rax,-0x38(%rbp)
   143c0:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   143c4:	48 8b 00             	mov    (%rax),%rax
   143c7:	48 8b 40 18          	mov    0x18(%rax),%rax
   143cb:	48 89 45 d0          	mov    %rax,-0x30(%rbp)
   143cf:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   143d3:	48 8b 00             	mov    (%rax),%rax
   143d6:	48 8b 40 20          	mov    0x20(%rax),%rax
   143da:	48 89 45 d8          	mov    %rax,-0x28(%rbp)
   143de:	48 8b 45 10          	mov    0x10(%rbp),%rax
   143e2:	48 8b 48 28          	mov    0x28(%rax),%rcx
   143e6:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   143ea:	8b 00                	mov    (%rax),%eax
   143ec:	48 98                	cltq   
   143ee:	48 c1 e0 06          	shl    $0x6,%rax
   143f2:	48 89 c6             	mov    %rax,%rsi
   143f5:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   143f9:	8b 00                	mov    (%rax),%eax
   143fb:	48 63 f8             	movslq %eax,%rdi
   143fe:	48 89 f0             	mov    %rsi,%rax
   14401:	ba 00 00 00 00       	mov    $0x0,%edx
   14406:	48 f7 f7             	div    %rdi
   14409:	89 c7                	mov    %eax,%edi
   1440b:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   1440f:	48 8b 10             	mov    (%rax),%rdx
   14412:	8b 45 a0             	mov    -0x60(%rbp),%eax
   14415:	48 98                	cltq   
   14417:	c1 e0 06             	shl    $0x6,%eax
   1441a:	89 c6                	mov    %eax,%esi
   1441c:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   14420:	51                   	push   %rcx
   14421:	6a 00                	push   $0x0
   14423:	4c 8d 0d f6 93 00 00 	lea    0x93f6(%rip),%r9        # 1d820 <ompi_mpi_byte@@Base>
   1442a:	41 89 f8             	mov    %edi,%r8d
   1442d:	48 89 d1             	mov    %rdx,%rcx
   14430:	48 8d 15 e9 93 00 00 	lea    0x93e9(%rip),%rdx        # 1d820 <ompi_mpi_byte@@Base>
   14437:	48 89 c7             	mov    %rax,%rdi
   1443a:	e8 41 5e ff ff       	call   a280 <MPI_Gather@plt>
   1443f:	48 83 c4 10          	add    $0x10,%rsp
   14443:	48 8b 45 10          	mov    0x10(%rbp),%rax
   14447:	48 8b 48 28          	mov    0x28(%rax),%rcx
   1444b:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   1444f:	8b 10                	mov    (%rax),%edx
   14451:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   14455:	48 8b 00             	mov    (%rax),%rax
   14458:	8b 40 30             	mov    0x30(%rax),%eax
   1445b:	0f af c2             	imul   %edx,%eax
   1445e:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   14462:	8b 3a                	mov    (%rdx),%edi
   14464:	99                   	cltd   
   14465:	f7 ff                	idiv   %edi
   14467:	89 c7                	mov    %eax,%edi
   14469:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1446d:	8b 40 30             	mov    0x30(%rax),%eax
   14470:	0f af 45 a0          	imul   -0x60(%rbp),%eax
   14474:	89 c6                	mov    %eax,%esi
   14476:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1447a:	48 8b 40 08          	mov    0x8(%rax),%rax
   1447e:	48 8b 55 c0          	mov    -0x40(%rbp),%rdx
   14482:	51                   	push   %rcx
   14483:	6a 00                	push   $0x0
   14485:	4c 8d 0d 94 97 00 00 	lea    0x9794(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   1448c:	41 89 f8             	mov    %edi,%r8d
   1448f:	48 89 d1             	mov    %rdx,%rcx
   14492:	48 8d 15 87 97 00 00 	lea    0x9787(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   14499:	48 89 c7             	mov    %rax,%rdi
   1449c:	e8 df 5d ff ff       	call   a280 <MPI_Gather@plt>
   144a1:	48 83 c4 10          	add    $0x10,%rsp
   144a5:	48 8b 45 10          	mov    0x10(%rbp),%rax
   144a9:	48 8b 48 28          	mov    0x28(%rax),%rcx
   144ad:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   144b1:	8b 10                	mov    (%rax),%edx
   144b3:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   144b7:	48 8b 00             	mov    (%rax),%rax
   144ba:	8b 40 30             	mov    0x30(%rax),%eax
   144bd:	0f af c2             	imul   %edx,%eax
   144c0:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   144c4:	8b 3a                	mov    (%rdx),%edi
   144c6:	99                   	cltd   
   144c7:	f7 ff                	idiv   %edi
   144c9:	89 c7                	mov    %eax,%edi
   144cb:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   144cf:	8b 40 30             	mov    0x30(%rax),%eax
   144d2:	0f af 45 a0          	imul   -0x60(%rbp),%eax
   144d6:	89 c6                	mov    %eax,%esi
   144d8:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   144dc:	48 8b 40 10          	mov    0x10(%rax),%rax
   144e0:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
   144e4:	51                   	push   %rcx
   144e5:	6a 00                	push   $0x0
   144e7:	4c 8d 0d 32 97 00 00 	lea    0x9732(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   144ee:	41 89 f8             	mov    %edi,%r8d
   144f1:	48 89 d1             	mov    %rdx,%rcx
   144f4:	48 8d 15 25 97 00 00 	lea    0x9725(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   144fb:	48 89 c7             	mov    %rax,%rdi
   144fe:	e8 7d 5d ff ff       	call   a280 <MPI_Gather@plt>
   14503:	48 83 c4 10          	add    $0x10,%rsp
   14507:	48 8b 45 10          	mov    0x10(%rbp),%rax
   1450b:	48 8b 48 28          	mov    0x28(%rax),%rcx
   1450f:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14513:	8b 10                	mov    (%rax),%edx
   14515:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   14519:	48 8b 00             	mov    (%rax),%rax
   1451c:	8b 40 30             	mov    0x30(%rax),%eax
   1451f:	0f af c2             	imul   %edx,%eax
   14522:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   14526:	8b 3a                	mov    (%rdx),%edi
   14528:	99                   	cltd   
   14529:	f7 ff                	idiv   %edi
   1452b:	89 c7                	mov    %eax,%edi
   1452d:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   14531:	8b 40 30             	mov    0x30(%rax),%eax
   14534:	0f af 45 a0          	imul   -0x60(%rbp),%eax
   14538:	89 c6                	mov    %eax,%esi
   1453a:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1453e:	48 8b 40 18          	mov    0x18(%rax),%rax
   14542:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
   14546:	51                   	push   %rcx
   14547:	6a 00                	push   $0x0
   14549:	4c 8d 0d d0 96 00 00 	lea    0x96d0(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   14550:	41 89 f8             	mov    %edi,%r8d
   14553:	48 89 d1             	mov    %rdx,%rcx
   14556:	48 8d 15 c3 96 00 00 	lea    0x96c3(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   1455d:	48 89 c7             	mov    %rax,%rdi
   14560:	e8 1b 5d ff ff       	call   a280 <MPI_Gather@plt>
   14565:	48 83 c4 10          	add    $0x10,%rsp
   14569:	48 8b 45 10          	mov    0x10(%rbp),%rax
   1456d:	48 8b 48 28          	mov    0x28(%rax),%rcx
   14571:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14575:	8b 10                	mov    (%rax),%edx
   14577:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   1457b:	48 8b 00             	mov    (%rax),%rax
   1457e:	8b 40 34             	mov    0x34(%rax),%eax
   14581:	0f af c2             	imul   %edx,%eax
   14584:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   14588:	8b 3a                	mov    (%rdx),%edi
   1458a:	99                   	cltd   
   1458b:	f7 ff                	idiv   %edi
   1458d:	89 c7                	mov    %eax,%edi
   1458f:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   14593:	8b 40 34             	mov    0x34(%rax),%eax
   14596:	0f af 45 a0          	imul   -0x60(%rbp),%eax
   1459a:	89 c6                	mov    %eax,%esi
   1459c:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   145a0:	48 8b 40 20          	mov    0x20(%rax),%rax
   145a4:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
   145a8:	51                   	push   %rcx
   145a9:	6a 00                	push   $0x0
   145ab:	4c 8d 0d 6e 96 00 00 	lea    0x966e(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   145b2:	41 89 f8             	mov    %edi,%r8d
   145b5:	48 89 d1             	mov    %rdx,%rcx
   145b8:	48 8d 15 61 96 00 00 	lea    0x9661(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   145bf:	48 89 c7             	mov    %rax,%rdi
   145c2:	e8 b9 5c ff ff       	call   a280 <MPI_Gather@plt>
   145c7:	48 83 c4 10          	add    $0x10,%rsp
   145cb:	e9 52 01 00 00       	jmp    14722 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x47a>
   145d0:	48 8b 45 10          	mov    0x10(%rbp),%rax
   145d4:	48 8b 50 28          	mov    0x28(%rax),%rdx
   145d8:	8b 45 a0             	mov    -0x60(%rbp),%eax
   145db:	48 98                	cltq   
   145dd:	c1 e0 06             	shl    $0x6,%eax
   145e0:	89 c6                	mov    %eax,%esi
   145e2:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   145e6:	52                   	push   %rdx
   145e7:	6a 00                	push   $0x0
   145e9:	4c 8d 0d 30 92 00 00 	lea    0x9230(%rip),%r9        # 1d820 <ompi_mpi_byte@@Base>
   145f0:	41 b8 ff ff ff ff    	mov    $0xffffffff,%r8d
   145f6:	b9 00 00 00 00       	mov    $0x0,%ecx
   145fb:	48 8d 15 1e 92 00 00 	lea    0x921e(%rip),%rdx        # 1d820 <ompi_mpi_byte@@Base>
   14602:	48 89 c7             	mov    %rax,%rdi
   14605:	e8 76 5c ff ff       	call   a280 <MPI_Gather@plt>
   1460a:	48 83 c4 10          	add    $0x10,%rsp
   1460e:	48 8b 45 10          	mov    0x10(%rbp),%rax
   14612:	48 8b 50 28          	mov    0x28(%rax),%rdx
   14616:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1461a:	8b 40 30             	mov    0x30(%rax),%eax
   1461d:	0f af 45 a0          	imul   -0x60(%rbp),%eax
   14621:	89 c6                	mov    %eax,%esi
   14623:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   14627:	48 8b 40 08          	mov    0x8(%rax),%rax
   1462b:	52                   	push   %rdx
   1462c:	6a 00                	push   $0x0
   1462e:	4c 8d 0d eb 95 00 00 	lea    0x95eb(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   14635:	41 b8 ff ff ff ff    	mov    $0xffffffff,%r8d
   1463b:	b9 00 00 00 00       	mov    $0x0,%ecx
   14640:	48 8d 15 d9 95 00 00 	lea    0x95d9(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   14647:	48 89 c7             	mov    %rax,%rdi
   1464a:	e8 31 5c ff ff       	call   a280 <MPI_Gather@plt>
   1464f:	48 83 c4 10          	add    $0x10,%rsp
   14653:	48 8b 45 10          	mov    0x10(%rbp),%rax
   14657:	48 8b 50 28          	mov    0x28(%rax),%rdx
   1465b:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1465f:	8b 40 30             	mov    0x30(%rax),%eax
   14662:	0f af 45 a0          	imul   -0x60(%rbp),%eax
   14666:	89 c6                	mov    %eax,%esi
   14668:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   1466c:	48 8b 40 10          	mov    0x10(%rax),%rax
   14670:	52                   	push   %rdx
   14671:	6a 00                	push   $0x0
   14673:	4c 8d 0d a6 95 00 00 	lea    0x95a6(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   1467a:	41 b8 ff ff ff ff    	mov    $0xffffffff,%r8d
   14680:	b9 00 00 00 00       	mov    $0x0,%ecx
   14685:	48 8d 15 94 95 00 00 	lea    0x9594(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   1468c:	48 89 c7             	mov    %rax,%rdi
   1468f:	e8 ec 5b ff ff       	call   a280 <MPI_Gather@plt>
   14694:	48 83 c4 10          	add    $0x10,%rsp
   14698:	48 8b 45 10          	mov    0x10(%rbp),%rax
   1469c:	48 8b 50 28          	mov    0x28(%rax),%rdx
   146a0:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   146a4:	8b 40 30             	mov    0x30(%rax),%eax
   146a7:	0f af 45 a0          	imul   -0x60(%rbp),%eax
   146ab:	89 c6                	mov    %eax,%esi
   146ad:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   146b1:	48 8b 40 18          	mov    0x18(%rax),%rax
   146b5:	52                   	push   %rdx
   146b6:	6a 00                	push   $0x0
   146b8:	4c 8d 0d 61 95 00 00 	lea    0x9561(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   146bf:	41 b8 ff ff ff ff    	mov    $0xffffffff,%r8d
   146c5:	b9 00 00 00 00       	mov    $0x0,%ecx
   146ca:	48 8d 15 4f 95 00 00 	lea    0x954f(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   146d1:	48 89 c7             	mov    %rax,%rdi
   146d4:	e8 a7 5b ff ff       	call   a280 <MPI_Gather@plt>
   146d9:	48 83 c4 10          	add    $0x10,%rsp
   146dd:	48 8b 45 10          	mov    0x10(%rbp),%rax
   146e1:	48 8b 50 28          	mov    0x28(%rax),%rdx
   146e5:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   146e9:	8b 40 34             	mov    0x34(%rax),%eax
   146ec:	0f af 45 a0          	imul   -0x60(%rbp),%eax
   146f0:	89 c6                	mov    %eax,%esi
   146f2:	48 8b 45 98          	mov    -0x68(%rbp),%rax
   146f6:	48 8b 40 20          	mov    0x20(%rax),%rax
   146fa:	52                   	push   %rdx
   146fb:	6a 00                	push   $0x0
   146fd:	4c 8d 0d 1c 95 00 00 	lea    0x951c(%rip),%r9        # 1dc20 <ompi_mpi_double@@Base>
   14704:	41 b8 ff ff ff ff    	mov    $0xffffffff,%r8d
   1470a:	b9 00 00 00 00       	mov    $0x0,%ecx
   1470f:	48 8d 15 0a 95 00 00 	lea    0x950a(%rip),%rdx        # 1dc20 <ompi_mpi_double@@Base>
   14716:	48 89 c7             	mov    %rax,%rdi
   14719:	e8 62 5b ff ff       	call   a280 <MPI_Gather@plt>
   1471e:	48 83 c4 10          	add    $0x10,%rsp
   14722:	48 8d 05 f7 92 00 00 	lea    0x92f7(%rip),%rax        # 1da20 <ompi_mpi_comm_world@@Base>
   14729:	48 89 c7             	mov    %rax,%rdi
   1472c:	e8 8f 5a ff ff       	call   a1c0 <MPI_Barrier@plt>
   14731:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   14735:	8b 00                	mov    (%rax),%eax
   14737:	85 c0                	test   %eax,%eax
   14739:	0f 85 06 01 00 00    	jne    14845 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x59d>
   1473f:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   14743:	8b 00                	mov    (%rax),%eax
   14745:	85 c0                	test   %eax,%eax
   14747:	0f 85 f8 00 00 00    	jne    14845 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x59d>
   1474d:	c7 45 bc 00 00 00 00 	movl   $0x0,-0x44(%rbp)
   14754:	e9 dc 00 00 00       	jmp    14835 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x58d>
   14759:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   1475d:	48 8b 00             	mov    (%rax),%rax
   14760:	8b 40 30             	mov    0x30(%rax),%eax
   14763:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   14767:	48 98                	cltq   
   14769:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   14770:	00 
   14771:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   14775:	48 8b 10             	mov    (%rax),%rdx
   14778:	8b 45 bc             	mov    -0x44(%rbp),%eax
   1477b:	48 98                	cltq   
   1477d:	48 c1 e0 06          	shl    $0x6,%rax
   14781:	48 01 d0             	add    %rdx,%rax
   14784:	48 8b 55 c0          	mov    -0x40(%rbp),%rdx
   14788:	48 01 ca             	add    %rcx,%rdx
   1478b:	48 89 50 08          	mov    %rdx,0x8(%rax)
   1478f:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   14793:	48 8b 00             	mov    (%rax),%rax
   14796:	8b 40 30             	mov    0x30(%rax),%eax
   14799:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   1479d:	48 98                	cltq   
   1479f:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   147a6:	00 
   147a7:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   147ab:	48 8b 10             	mov    (%rax),%rdx
   147ae:	8b 45 bc             	mov    -0x44(%rbp),%eax
   147b1:	48 98                	cltq   
   147b3:	48 c1 e0 06          	shl    $0x6,%rax
   147b7:	48 01 d0             	add    %rdx,%rax
   147ba:	48 8b 55 c8          	mov    -0x38(%rbp),%rdx
   147be:	48 01 ca             	add    %rcx,%rdx
   147c1:	48 89 50 10          	mov    %rdx,0x10(%rax)
   147c5:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   147c9:	48 8b 00             	mov    (%rax),%rax
   147cc:	8b 40 30             	mov    0x30(%rax),%eax
   147cf:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   147d3:	48 98                	cltq   
   147d5:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   147dc:	00 
   147dd:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   147e1:	48 8b 10             	mov    (%rax),%rdx
   147e4:	8b 45 bc             	mov    -0x44(%rbp),%eax
   147e7:	48 98                	cltq   
   147e9:	48 c1 e0 06          	shl    $0x6,%rax
   147ed:	48 01 d0             	add    %rdx,%rax
   147f0:	48 8b 55 d0          	mov    -0x30(%rbp),%rdx
   147f4:	48 01 ca             	add    %rcx,%rdx
   147f7:	48 89 50 18          	mov    %rdx,0x18(%rax)
   147fb:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   147ff:	48 8b 00             	mov    (%rax),%rax
   14802:	8b 40 34             	mov    0x34(%rax),%eax
   14805:	0f af 45 bc          	imul   -0x44(%rbp),%eax
   14809:	48 98                	cltq   
   1480b:	48 8d 0c c5 00 00 00 	lea    0x0(,%rax,8),%rcx
   14812:	00 
   14813:	48 8b 45 90          	mov    -0x70(%rbp),%rax
   14817:	48 8b 10             	mov    (%rax),%rdx
   1481a:	8b 45 bc             	mov    -0x44(%rbp),%eax
   1481d:	48 98                	cltq   
   1481f:	48 c1 e0 06          	shl    $0x6,%rax
   14823:	48 01 d0             	add    %rdx,%rax
   14826:	48 8b 55 d8          	mov    -0x28(%rbp),%rdx
   1482a:	48 01 ca             	add    %rcx,%rdx
   1482d:	48 89 50 20          	mov    %rdx,0x20(%rax)
   14831:	83 45 bc 01          	addl   $0x1,-0x44(%rbp)
   14835:	48 8b 45 88          	mov    -0x78(%rbp),%rax
   14839:	8b 00                	mov    (%rax),%eax
   1483b:	39 45 bc             	cmp    %eax,-0x44(%rbp)
   1483e:	0f 8c 15 ff ff ff    	jl     14759 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t+0x4b1>
   14844:	90                   	nop
   14845:	90                   	nop
   14846:	c9                   	leave  
   14847:	c3                   	ret    

0000000000014848 <_Z6gatherP17chunk_particles_tS0_S0_P10topology_t>:
   14848:	f3 0f 1e fa          	endbr64 
   1484c:	55                   	push   %rbp
   1484d:	48 89 e5             	mov    %rsp,%rbp
   14850:	48 83 ec 50          	sub    $0x50,%rsp
   14854:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
   14858:	48 89 75 c0          	mov    %rsi,-0x40(%rbp)
   1485c:	48 89 55 b8          	mov    %rdx,-0x48(%rbp)
   14860:	48 89 4d b0          	mov    %rcx,-0x50(%rbp)
   14864:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   1486b:	00 00 
   1486d:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   14871:	31 c0                	xor    %eax,%eax
   14873:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14877:	48 8b 00             	mov    (%rax),%rax
   1487a:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   1487e:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14882:	8b 40 08             	mov    0x8(%rax),%eax
   14885:	89 45 d8             	mov    %eax,-0x28(%rbp)
   14888:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   1488c:	48 8b 00             	mov    (%rax),%rax
   1488f:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   14893:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   14897:	8b 40 08             	mov    0x8(%rax),%eax
   1489a:	89 45 dc             	mov    %eax,-0x24(%rbp)
   1489d:	48 83 7d b8 00       	cmpq   $0x0,-0x48(%rbp)
   148a2:	75 11                	jne    148b5 <_Z6gatherP17chunk_particles_tS0_S0_P10topology_t+0x6d>
   148a4:	48 c7 45 e0 00 00 00 	movq   $0x0,-0x20(%rbp)
   148ab:	00 
   148ac:	c7 45 d4 ff ff ff ff 	movl   $0xffffffff,-0x2c(%rbp)
   148b3:	eb 45                	jmp    148fa <_Z6gatherP17chunk_particles_tS0_S0_P10topology_t+0xb2>
   148b5:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   148b9:	48 8b 00             	mov    (%rax),%rax
   148bc:	48 89 45 e0          	mov    %rax,-0x20(%rbp)
   148c0:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   148c4:	8b 40 08             	mov    0x8(%rax),%eax
   148c7:	89 45 d4             	mov    %eax,-0x2c(%rbp)
   148ca:	8b 55 d4             	mov    -0x2c(%rbp),%edx
   148cd:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   148d1:	41 89 d0             	mov    %edx,%r8d
   148d4:	48 89 c1             	mov    %rax,%rcx
   148d7:	ba c8 02 00 00       	mov    $0x2c8,%edx
   148dc:	48 8d 05 9d 23 00 00 	lea    0x239d(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   148e3:	48 89 c6             	mov    %rax,%rsi
   148e6:	48 8d 05 49 26 00 00 	lea    0x2649(%rip),%rax        # 16f36 <_ZN3MPIL11LOCK_SHAREDE+0x2ee>
   148ed:	48 89 c7             	mov    %rax,%rdi
   148f0:	b8 00 00 00 00       	mov    $0x0,%eax
   148f5:	e8 c6 54 ff ff       	call   9dc0 <printf@plt>
   148fa:	4c 8d 45 d4          	lea    -0x2c(%rbp),%r8
   148fe:	48 8d 7d e0          	lea    -0x20(%rbp),%rdi
   14902:	8b 4d dc             	mov    -0x24(%rbp),%ecx
   14905:	48 8b 55 f0          	mov    -0x10(%rbp),%rdx
   14909:	8b 75 d8             	mov    -0x28(%rbp),%esi
   1490c:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   14910:	48 83 ec 08          	sub    $0x8,%rsp
   14914:	ff 75 b0             	push   -0x50(%rbp)
   14917:	4d 89 c1             	mov    %r8,%r9
   1491a:	49 89 f8             	mov    %rdi,%r8
   1491d:	48 89 c7             	mov    %rax,%rdi
   14920:	e8 83 f9 ff ff       	call   142a8 <_Z6gatherP10particle_tiS0_iPS0_PiP10topology_t>
   14925:	48 83 c4 10          	add    $0x10,%rsp
   14929:	48 83 7d b8 00       	cmpq   $0x0,-0x48(%rbp)
   1492e:	74 45                	je     14975 <_Z6gatherP17chunk_particles_tS0_S0_P10topology_t+0x12d>
   14930:	8b 55 d4             	mov    -0x2c(%rbp),%edx
   14933:	48 8b 45 e0          	mov    -0x20(%rbp),%rax
   14937:	41 89 d0             	mov    %edx,%r8d
   1493a:	48 89 c1             	mov    %rax,%rcx
   1493d:	ba ce 02 00 00       	mov    $0x2ce,%edx
   14942:	48 8d 05 37 23 00 00 	lea    0x2337(%rip),%rax        # 16c80 <_ZN3MPIL11LOCK_SHAREDE+0x38>
   14949:	48 89 c6             	mov    %rax,%rsi
   1494c:	48 8d 05 e3 25 00 00 	lea    0x25e3(%rip),%rax        # 16f36 <_ZN3MPIL11LOCK_SHAREDE+0x2ee>
   14953:	48 89 c7             	mov    %rax,%rdi
   14956:	b8 00 00 00 00       	mov    $0x0,%eax
   1495b:	e8 60 54 ff ff       	call   9dc0 <printf@plt>
   14960:	48 8b 55 e0          	mov    -0x20(%rbp),%rdx
   14964:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14968:	48 89 10             	mov    %rdx,(%rax)
   1496b:	8b 55 d4             	mov    -0x2c(%rbp),%edx
   1496e:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14972:	89 50 08             	mov    %edx,0x8(%rax)
   14975:	90                   	nop
   14976:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1497a:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   14981:	00 00 
   14983:	74 05                	je     1498a <_Z6gatherP17chunk_particles_tS0_S0_P10topology_t+0x142>
   14985:	e8 46 5b ff ff       	call   a4d0 <__stack_chk_fail@plt>
   1498a:	c9                   	leave  
   1498b:	c3                   	ret    

000000000001498c <_Z6updateP10particle_tiS0_iP10topology_td>:
   1498c:	f3 0f 1e fa          	endbr64 
   14990:	55                   	push   %rbp
   14991:	48 89 e5             	mov    %rsp,%rbp
   14994:	48 83 ec 60          	sub    $0x60,%rsp
   14998:	48 89 7d c8          	mov    %rdi,-0x38(%rbp)
   1499c:	89 75 c4             	mov    %esi,-0x3c(%rbp)
   1499f:	48 89 55 b8          	mov    %rdx,-0x48(%rbp)
   149a3:	89 4d c0             	mov    %ecx,-0x40(%rbp)
   149a6:	4c 89 45 b0          	mov    %r8,-0x50(%rbp)
   149aa:	f2 0f 11 45 a8       	movsd  %xmm0,-0x58(%rbp)
   149af:	c7 45 d8 00 00 00 00 	movl   $0x0,-0x28(%rbp)
   149b6:	eb 48                	jmp    14a00 <_Z6updateP10particle_tiS0_iP10topology_td+0x74>
   149b8:	c7 45 dc 00 00 00 00 	movl   $0x0,-0x24(%rbp)
   149bf:	eb 2f                	jmp    149f0 <_Z6updateP10particle_tiS0_iP10topology_td+0x64>
   149c1:	8b 45 d8             	mov    -0x28(%rbp),%eax
   149c4:	48 98                	cltq   
   149c6:	48 c1 e0 06          	shl    $0x6,%rax
   149ca:	48 89 c2             	mov    %rax,%rdx
   149cd:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   149d1:	48 01 d0             	add    %rdx,%rax
   149d4:	48 8b 50 18          	mov    0x18(%rax),%rdx
   149d8:	8b 45 dc             	mov    -0x24(%rbp),%eax
   149db:	48 98                	cltq   
   149dd:	48 c1 e0 03          	shl    $0x3,%rax
   149e1:	48 01 d0             	add    %rdx,%rax
   149e4:	66 0f ef c0          	pxor   %xmm0,%xmm0
   149e8:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   149ec:	83 45 dc 01          	addl   $0x1,-0x24(%rbp)
   149f0:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   149f4:	8b 40 30             	mov    0x30(%rax),%eax
   149f7:	39 45 dc             	cmp    %eax,-0x24(%rbp)
   149fa:	7c c5                	jl     149c1 <_Z6updateP10particle_tiS0_iP10topology_td+0x35>
   149fc:	83 45 d8 01          	addl   $0x1,-0x28(%rbp)
   14a00:	8b 45 d8             	mov    -0x28(%rbp),%eax
   14a03:	3b 45 c4             	cmp    -0x3c(%rbp),%eax
   14a06:	7c b0                	jl     149b8 <_Z6updateP10particle_tiS0_iP10topology_td+0x2c>
   14a08:	c7 45 e0 00 00 00 00 	movl   $0x0,-0x20(%rbp)
   14a0f:	eb 48                	jmp    14a59 <_Z6updateP10particle_tiS0_iP10topology_td+0xcd>
   14a11:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
   14a18:	eb 2f                	jmp    14a49 <_Z6updateP10particle_tiS0_iP10topology_td+0xbd>
   14a1a:	8b 45 e0             	mov    -0x20(%rbp),%eax
   14a1d:	48 98                	cltq   
   14a1f:	48 c1 e0 06          	shl    $0x6,%rax
   14a23:	48 89 c2             	mov    %rax,%rdx
   14a26:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14a2a:	48 01 d0             	add    %rdx,%rax
   14a2d:	48 8b 50 18          	mov    0x18(%rax),%rdx
   14a31:	8b 45 e4             	mov    -0x1c(%rbp),%eax
   14a34:	48 98                	cltq   
   14a36:	48 c1 e0 03          	shl    $0x3,%rax
   14a3a:	48 01 d0             	add    %rdx,%rax
   14a3d:	66 0f ef c0          	pxor   %xmm0,%xmm0
   14a41:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   14a45:	83 45 e4 01          	addl   $0x1,-0x1c(%rbp)
   14a49:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14a4d:	8b 40 30             	mov    0x30(%rax),%eax
   14a50:	39 45 e4             	cmp    %eax,-0x1c(%rbp)
   14a53:	7c c5                	jl     14a1a <_Z6updateP10particle_tiS0_iP10topology_td+0x8e>
   14a55:	83 45 e0 01          	addl   $0x1,-0x20(%rbp)
   14a59:	8b 45 e0             	mov    -0x20(%rbp),%eax
   14a5c:	3b 45 c0             	cmp    -0x40(%rbp),%eax
   14a5f:	7c b0                	jl     14a11 <_Z6updateP10particle_tiS0_iP10topology_td+0x85>
   14a61:	c7 45 e8 00 00 00 00 	movl   $0x0,-0x18(%rbp)
   14a68:	eb 4a                	jmp    14ab4 <_Z6updateP10particle_tiS0_iP10topology_td+0x128>
   14a6a:	c7 45 ec 00 00 00 00 	movl   $0x0,-0x14(%rbp)
   14a71:	eb 35                	jmp    14aa8 <_Z6updateP10particle_tiS0_iP10topology_td+0x11c>
   14a73:	8b 45 ec             	mov    -0x14(%rbp),%eax
   14a76:	48 98                	cltq   
   14a78:	48 c1 e0 06          	shl    $0x6,%rax
   14a7c:	48 89 c2             	mov    %rax,%rdx
   14a7f:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14a83:	48 01 c2             	add    %rax,%rdx
   14a86:	8b 45 e8             	mov    -0x18(%rbp),%eax
   14a89:	48 98                	cltq   
   14a8b:	48 c1 e0 06          	shl    $0x6,%rax
   14a8f:	48 89 c1             	mov    %rax,%rcx
   14a92:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14a96:	48 01 c8             	add    %rcx,%rax
   14a99:	48 89 d6             	mov    %rdx,%rsi
   14a9c:	48 89 c7             	mov    %rax,%rdi
   14a9f:	e8 7d c0 ff ff       	call   10b21 <_Z18updateAccelerationP10particle_tS0_>
   14aa4:	83 45 ec 01          	addl   $0x1,-0x14(%rbp)
   14aa8:	8b 45 ec             	mov    -0x14(%rbp),%eax
   14aab:	3b 45 c0             	cmp    -0x40(%rbp),%eax
   14aae:	7c c3                	jl     14a73 <_Z6updateP10particle_tiS0_iP10topology_td+0xe7>
   14ab0:	83 45 e8 01          	addl   $0x1,-0x18(%rbp)
   14ab4:	8b 45 e8             	mov    -0x18(%rbp),%eax
   14ab7:	3b 45 c4             	cmp    -0x3c(%rbp),%eax
   14aba:	7c ae                	jl     14a6a <_Z6updateP10particle_tiS0_iP10topology_td+0xde>
   14abc:	48 8b 7d b0          	mov    -0x50(%rbp),%rdi
   14ac0:	8b 4d c0             	mov    -0x40(%rbp),%ecx
   14ac3:	48 8b 55 b8          	mov    -0x48(%rbp),%rdx
   14ac7:	8b 75 c4             	mov    -0x3c(%rbp),%esi
   14aca:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14ace:	49 89 f8             	mov    %rdi,%r8
   14ad1:	48 89 c7             	mov    %rax,%rdi
   14ad4:	e8 a2 f5 ff ff       	call   1407b <_Z6reduceP10particle_tiS0_iP10topology_t>
   14ad9:	c7 45 f0 00 00 00 00 	movl   $0x0,-0x10(%rbp)
   14ae0:	e9 c6 00 00 00       	jmp    14bab <_Z6updateP10particle_tiS0_iP10topology_td+0x21f>
   14ae5:	c7 45 f4 00 00 00 00 	movl   $0x0,-0xc(%rbp)
   14aec:	e9 a6 00 00 00       	jmp    14b97 <_Z6updateP10particle_tiS0_iP10topology_td+0x20b>
   14af1:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14af5:	48 8b 50 10          	mov    0x10(%rax),%rdx
   14af9:	8b 45 f4             	mov    -0xc(%rbp),%eax
   14afc:	48 98                	cltq   
   14afe:	48 c1 e0 03          	shl    $0x3,%rax
   14b02:	48 01 d0             	add    %rdx,%rax
   14b05:	f2 0f 10 08          	movsd  (%rax),%xmm1
   14b09:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14b0d:	48 8b 50 18          	mov    0x18(%rax),%rdx
   14b11:	8b 45 f4             	mov    -0xc(%rbp),%eax
   14b14:	48 98                	cltq   
   14b16:	48 c1 e0 03          	shl    $0x3,%rax
   14b1a:	48 01 d0             	add    %rdx,%rax
   14b1d:	f2 0f 10 00          	movsd  (%rax),%xmm0
   14b21:	f2 0f 59 45 a8       	mulsd  -0x58(%rbp),%xmm0
   14b26:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14b2a:	48 8b 50 10          	mov    0x10(%rax),%rdx
   14b2e:	8b 45 f4             	mov    -0xc(%rbp),%eax
   14b31:	48 98                	cltq   
   14b33:	48 c1 e0 03          	shl    $0x3,%rax
   14b37:	48 01 d0             	add    %rdx,%rax
   14b3a:	f2 0f 58 c1          	addsd  %xmm1,%xmm0
   14b3e:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   14b42:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14b46:	48 8b 50 08          	mov    0x8(%rax),%rdx
   14b4a:	8b 45 f4             	mov    -0xc(%rbp),%eax
   14b4d:	48 98                	cltq   
   14b4f:	48 c1 e0 03          	shl    $0x3,%rax
   14b53:	48 01 d0             	add    %rdx,%rax
   14b56:	f2 0f 10 08          	movsd  (%rax),%xmm1
   14b5a:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14b5e:	48 8b 50 10          	mov    0x10(%rax),%rdx
   14b62:	8b 45 f4             	mov    -0xc(%rbp),%eax
   14b65:	48 98                	cltq   
   14b67:	48 c1 e0 03          	shl    $0x3,%rax
   14b6b:	48 01 d0             	add    %rdx,%rax
   14b6e:	f2 0f 10 00          	movsd  (%rax),%xmm0
   14b72:	f2 0f 59 45 a8       	mulsd  -0x58(%rbp),%xmm0
   14b77:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14b7b:	48 8b 50 08          	mov    0x8(%rax),%rdx
   14b7f:	8b 45 f4             	mov    -0xc(%rbp),%eax
   14b82:	48 98                	cltq   
   14b84:	48 c1 e0 03          	shl    $0x3,%rax
   14b88:	48 01 d0             	add    %rdx,%rax
   14b8b:	f2 0f 58 c1          	addsd  %xmm1,%xmm0
   14b8f:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   14b93:	83 45 f4 01          	addl   $0x1,-0xc(%rbp)
   14b97:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14b9b:	8b 40 30             	mov    0x30(%rax),%eax
   14b9e:	39 45 f4             	cmp    %eax,-0xc(%rbp)
   14ba1:	0f 8c 4a ff ff ff    	jl     14af1 <_Z6updateP10particle_tiS0_iP10topology_td+0x165>
   14ba7:	83 45 f0 01          	addl   $0x1,-0x10(%rbp)
   14bab:	8b 45 f0             	mov    -0x10(%rbp),%eax
   14bae:	3b 45 c4             	cmp    -0x3c(%rbp),%eax
   14bb1:	0f 8c 2e ff ff ff    	jl     14ae5 <_Z6updateP10particle_tiS0_iP10topology_td+0x159>
   14bb7:	c7 45 f8 00 00 00 00 	movl   $0x0,-0x8(%rbp)
   14bbe:	e9 c6 00 00 00       	jmp    14c89 <_Z6updateP10particle_tiS0_iP10topology_td+0x2fd>
   14bc3:	c7 45 fc 00 00 00 00 	movl   $0x0,-0x4(%rbp)
   14bca:	e9 a6 00 00 00       	jmp    14c75 <_Z6updateP10particle_tiS0_iP10topology_td+0x2e9>
   14bcf:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14bd3:	48 8b 50 10          	mov    0x10(%rax),%rdx
   14bd7:	8b 45 fc             	mov    -0x4(%rbp),%eax
   14bda:	48 98                	cltq   
   14bdc:	48 c1 e0 03          	shl    $0x3,%rax
   14be0:	48 01 d0             	add    %rdx,%rax
   14be3:	f2 0f 10 08          	movsd  (%rax),%xmm1
   14be7:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14beb:	48 8b 50 18          	mov    0x18(%rax),%rdx
   14bef:	8b 45 fc             	mov    -0x4(%rbp),%eax
   14bf2:	48 98                	cltq   
   14bf4:	48 c1 e0 03          	shl    $0x3,%rax
   14bf8:	48 01 d0             	add    %rdx,%rax
   14bfb:	f2 0f 10 00          	movsd  (%rax),%xmm0
   14bff:	f2 0f 59 45 a8       	mulsd  -0x58(%rbp),%xmm0
   14c04:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14c08:	48 8b 50 10          	mov    0x10(%rax),%rdx
   14c0c:	8b 45 fc             	mov    -0x4(%rbp),%eax
   14c0f:	48 98                	cltq   
   14c11:	48 c1 e0 03          	shl    $0x3,%rax
   14c15:	48 01 d0             	add    %rdx,%rax
   14c18:	f2 0f 58 c1          	addsd  %xmm1,%xmm0
   14c1c:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   14c20:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14c24:	48 8b 50 08          	mov    0x8(%rax),%rdx
   14c28:	8b 45 fc             	mov    -0x4(%rbp),%eax
   14c2b:	48 98                	cltq   
   14c2d:	48 c1 e0 03          	shl    $0x3,%rax
   14c31:	48 01 d0             	add    %rdx,%rax
   14c34:	f2 0f 10 08          	movsd  (%rax),%xmm1
   14c38:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14c3c:	48 8b 50 10          	mov    0x10(%rax),%rdx
   14c40:	8b 45 fc             	mov    -0x4(%rbp),%eax
   14c43:	48 98                	cltq   
   14c45:	48 c1 e0 03          	shl    $0x3,%rax
   14c49:	48 01 d0             	add    %rdx,%rax
   14c4c:	f2 0f 10 00          	movsd  (%rax),%xmm0
   14c50:	f2 0f 59 45 a8       	mulsd  -0x58(%rbp),%xmm0
   14c55:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14c59:	48 8b 50 08          	mov    0x8(%rax),%rdx
   14c5d:	8b 45 fc             	mov    -0x4(%rbp),%eax
   14c60:	48 98                	cltq   
   14c62:	48 c1 e0 03          	shl    $0x3,%rax
   14c66:	48 01 d0             	add    %rdx,%rax
   14c69:	f2 0f 58 c1          	addsd  %xmm1,%xmm0
   14c6d:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   14c71:	83 45 fc 01          	addl   $0x1,-0x4(%rbp)
   14c75:	48 8b 45 b8          	mov    -0x48(%rbp),%rax
   14c79:	8b 40 30             	mov    0x30(%rax),%eax
   14c7c:	39 45 fc             	cmp    %eax,-0x4(%rbp)
   14c7f:	0f 8c 4a ff ff ff    	jl     14bcf <_Z6updateP10particle_tiS0_iP10topology_td+0x243>
   14c85:	83 45 f8 01          	addl   $0x1,-0x8(%rbp)
   14c89:	8b 45 f8             	mov    -0x8(%rbp),%eax
   14c8c:	3b 45 c0             	cmp    -0x40(%rbp),%eax
   14c8f:	0f 8c 2e ff ff ff    	jl     14bc3 <_Z6updateP10particle_tiS0_iP10topology_td+0x237>
   14c95:	90                   	nop
   14c96:	90                   	nop
   14c97:	c9                   	leave  
   14c98:	c3                   	ret    

0000000000014c99 <_Z6updateP17chunk_particles_tS0_P10topology_td>:
   14c99:	f3 0f 1e fa          	endbr64 
   14c9d:	55                   	push   %rbp
   14c9e:	48 89 e5             	mov    %rsp,%rbp
   14ca1:	48 83 ec 20          	sub    $0x20,%rsp
   14ca5:	48 89 7d f8          	mov    %rdi,-0x8(%rbp)
   14ca9:	48 89 75 f0          	mov    %rsi,-0x10(%rbp)
   14cad:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
   14cb1:	f2 0f 11 45 e0       	movsd  %xmm0,-0x20(%rbp)
   14cb6:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   14cba:	8b 48 08             	mov    0x8(%rax),%ecx
   14cbd:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   14cc1:	48 8b 10             	mov    (%rax),%rdx
   14cc4:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14cc8:	8b 70 08             	mov    0x8(%rax),%esi
   14ccb:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14ccf:	48 8b 00             	mov    (%rax),%rax
   14cd2:	48 8b 7d e0          	mov    -0x20(%rbp),%rdi
   14cd6:	4c 8b 45 e8          	mov    -0x18(%rbp),%r8
   14cda:	66 48 0f 6e c7       	movq   %rdi,%xmm0
   14cdf:	48 89 c7             	mov    %rax,%rdi
   14ce2:	e8 a5 fc ff ff       	call   1498c <_Z6updateP10particle_tiS0_iP10topology_td>
   14ce7:	90                   	nop
   14ce8:	c9                   	leave  
   14ce9:	c3                   	ret    

0000000000014cea <_Z41__static_initialization_and_destruction_0ii>:
   14cea:	f3 0f 1e fa          	endbr64 
   14cee:	55                   	push   %rbp
   14cef:	48 89 e5             	mov    %rsp,%rbp
   14cf2:	48 83 ec 10          	sub    $0x10,%rsp
   14cf6:	89 7d fc             	mov    %edi,-0x4(%rbp)
   14cf9:	89 75 f8             	mov    %esi,-0x8(%rbp)
   14cfc:	83 7d fc 01          	cmpl   $0x1,-0x4(%rbp)
   14d00:	75 3b                	jne    14d3d <_Z41__static_initialization_and_destruction_0ii+0x53>
   14d02:	81 7d f8 ff ff 00 00 	cmpl   $0xffff,-0x8(%rbp)
   14d09:	75 32                	jne    14d3d <_Z41__static_initialization_and_destruction_0ii+0x53>
   14d0b:	48 8d 05 19 95 00 00 	lea    0x9519(%rip),%rax        # 1e22b <_ZStL8__ioinit>
   14d12:	48 89 c7             	mov    %rax,%rdi
   14d15:	e8 96 5a ff ff       	call   a7b0 <_ZNSt8ios_base4InitC1Ev@plt>
   14d1a:	48 8d 05 e7 82 00 00 	lea    0x82e7(%rip),%rax        # 1d008 <__dso_handle>
   14d21:	48 89 c2             	mov    %rax,%rdx
   14d24:	48 8d 05 00 95 00 00 	lea    0x9500(%rip),%rax        # 1e22b <_ZStL8__ioinit>
   14d2b:	48 89 c6             	mov    %rax,%rsi
   14d2e:	48 8b 05 c3 82 00 00 	mov    0x82c3(%rip),%rax        # 1cff8 <_ZNSt8ios_base4InitD1Ev@GLIBCXX_3.4>
   14d35:	48 89 c7             	mov    %rax,%rdi
   14d38:	e8 e3 55 ff ff       	call   a320 <__cxa_atexit@plt>
   14d3d:	90                   	nop
   14d3e:	c9                   	leave  
   14d3f:	c3                   	ret    

0000000000014d40 <_GLOBAL__sub_I__Z18setupCommunicatorsP10topology_tii>:
   14d40:	f3 0f 1e fa          	endbr64 
   14d44:	55                   	push   %rbp
   14d45:	48 89 e5             	mov    %rsp,%rbp
   14d48:	be ff ff 00 00       	mov    $0xffff,%esi
   14d4d:	bf 01 00 00 00       	mov    $0x1,%edi
   14d52:	e8 93 ff ff ff       	call   14cea <_Z41__static_initialization_and_destruction_0ii>
   14d57:	5d                   	pop    %rbp
   14d58:	c3                   	ret    

0000000000014d59 <_Z16sizeAfterPaddingii>:
   14d59:	f3 0f 1e fa          	endbr64 
   14d5d:	55                   	push   %rbp
   14d5e:	48 89 e5             	mov    %rsp,%rbp
   14d61:	89 7d fc             	mov    %edi,-0x4(%rbp)
   14d64:	89 75 f8             	mov    %esi,-0x8(%rbp)
   14d67:	8b 55 fc             	mov    -0x4(%rbp),%edx
   14d6a:	8b 45 f8             	mov    -0x8(%rbp),%eax
   14d6d:	01 d0                	add    %edx,%eax
   14d6f:	83 e8 01             	sub    $0x1,%eax
   14d72:	99                   	cltd   
   14d73:	f7 7d f8             	idivl  -0x8(%rbp)
   14d76:	0f af 45 f8          	imul   -0x8(%rbp),%eax
   14d7a:	5d                   	pop    %rbp
   14d7b:	c3                   	ret    

0000000000014d7c <_Z16get_distance_vecPKdS0_iPPd>:
   14d7c:	f3 0f 1e fa          	endbr64 
   14d80:	55                   	push   %rbp
   14d81:	48 89 e5             	mov    %rsp,%rbp
   14d84:	48 83 ec 40          	sub    $0x40,%rsp
   14d88:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   14d8c:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
   14d90:	89 55 cc             	mov    %edx,-0x34(%rbp)
   14d93:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
   14d97:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   14d9b:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   14d9f:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   14da3:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   14da7:	8b 45 cc             	mov    -0x34(%rbp),%eax
   14daa:	48 98                	cltq   
   14dac:	48 c1 e0 03          	shl    $0x3,%rax
   14db0:	48 89 c7             	mov    %rax,%rdi
   14db3:	e8 68 58 ff ff       	call   a620 <malloc@plt>
   14db8:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   14dbc:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   14dc0:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   14dc4:	48 89 10             	mov    %rdx,(%rax)
   14dc7:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
   14dce:	eb 38                	jmp    14e08 <_Z16get_distance_vecPKdS0_iPPd+0x8c>
   14dd0:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   14dd4:	48 8d 50 08          	lea    0x8(%rax),%rdx
   14dd8:	48 89 55 f0          	mov    %rdx,-0x10(%rbp)
   14ddc:	f2 0f 10 00          	movsd  (%rax),%xmm0
   14de0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14de4:	48 8d 50 08          	lea    0x8(%rax),%rdx
   14de8:	48 89 55 f8          	mov    %rdx,-0x8(%rbp)
   14dec:	f2 0f 10 08          	movsd  (%rax),%xmm1
   14df0:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   14df4:	48 8d 50 08          	lea    0x8(%rax),%rdx
   14df8:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
   14dfc:	f2 0f 5c c1          	subsd  %xmm1,%xmm0
   14e00:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   14e04:	83 45 e4 01          	addl   $0x1,-0x1c(%rbp)
   14e08:	8b 45 e4             	mov    -0x1c(%rbp),%eax
   14e0b:	3b 45 cc             	cmp    -0x34(%rbp),%eax
   14e0e:	7c c0                	jl     14dd0 <_Z16get_distance_vecPKdS0_iPPd+0x54>
   14e10:	90                   	nop
   14e11:	90                   	nop
   14e12:	c9                   	leave  
   14e13:	c3                   	ret    

0000000000014e14 <_Z10get_lengthPKdi>:
   14e14:	f3 0f 1e fa          	endbr64 
   14e18:	55                   	push   %rbp
   14e19:	48 89 e5             	mov    %rsp,%rbp
   14e1c:	48 83 ec 30          	sub    $0x30,%rsp
   14e20:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   14e24:	89 75 d4             	mov    %esi,-0x2c(%rbp)
   14e27:	66 0f ef c0          	pxor   %xmm0,%xmm0
   14e2b:	f2 0f 11 45 f0       	movsd  %xmm0,-0x10(%rbp)
   14e30:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   14e34:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   14e38:	c7 45 ec 00 00 00 00 	movl   $0x0,-0x14(%rbp)
   14e3f:	eb 2b                	jmp    14e6c <_Z10get_lengthPKdi+0x58>
   14e41:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14e45:	f2 0f 10 08          	movsd  (%rax),%xmm1
   14e49:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14e4d:	f2 0f 10 00          	movsd  (%rax),%xmm0
   14e51:	f2 0f 59 c1          	mulsd  %xmm1,%xmm0
   14e55:	f2 0f 10 4d f0       	movsd  -0x10(%rbp),%xmm1
   14e5a:	f2 0f 58 c1          	addsd  %xmm1,%xmm0
   14e5e:	f2 0f 11 45 f0       	movsd  %xmm0,-0x10(%rbp)
   14e63:	48 83 45 f8 08       	addq   $0x8,-0x8(%rbp)
   14e68:	83 45 ec 01          	addl   $0x1,-0x14(%rbp)
   14e6c:	8b 45 ec             	mov    -0x14(%rbp),%eax
   14e6f:	3b 45 d4             	cmp    -0x2c(%rbp),%eax
   14e72:	7c cd                	jl     14e41 <_Z10get_lengthPKdi+0x2d>
   14e74:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   14e78:	66 48 0f 6e c0       	movq   %rax,%xmm0
   14e7d:	e8 de 53 ff ff       	call   a260 <sqrt@plt>
   14e82:	66 48 0f 7e c0       	movq   %xmm0,%rax
   14e87:	66 48 0f 6e c0       	movq   %rax,%xmm0
   14e8c:	c9                   	leave  
   14e8d:	c3                   	ret    

0000000000014e8e <_Z7l2_normPKdiPPd>:
   14e8e:	f3 0f 1e fa          	endbr64 
   14e92:	55                   	push   %rbp
   14e93:	48 89 e5             	mov    %rsp,%rbp
   14e96:	48 83 ec 40          	sub    $0x40,%rsp
   14e9a:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   14e9e:	89 75 d4             	mov    %esi,-0x2c(%rbp)
   14ea1:	48 89 55 c8          	mov    %rdx,-0x38(%rbp)
   14ea5:	8b 45 d4             	mov    -0x2c(%rbp),%eax
   14ea8:	48 98                	cltq   
   14eaa:	48 c1 e0 03          	shl    $0x3,%rax
   14eae:	48 89 c7             	mov    %rax,%rdi
   14eb1:	e8 6a 57 ff ff       	call   a620 <malloc@plt>
   14eb6:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   14eba:	48 8b 45 c8          	mov    -0x38(%rbp),%rax
   14ebe:	48 8b 55 e8          	mov    -0x18(%rbp),%rdx
   14ec2:	48 89 10             	mov    %rdx,(%rax)
   14ec5:	8b 55 d4             	mov    -0x2c(%rbp),%edx
   14ec8:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   14ecc:	89 d6                	mov    %edx,%esi
   14ece:	48 89 c7             	mov    %rax,%rdi
   14ed1:	e8 3e ff ff ff       	call   14e14 <_Z10get_lengthPKdi>
   14ed6:	66 48 0f 7e c0       	movq   %xmm0,%rax
   14edb:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   14edf:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   14ee3:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   14ee7:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
   14eee:	eb 21                	jmp    14f11 <_Z7l2_normPKdiPPd+0x83>
   14ef0:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14ef4:	f2 0f 10 00          	movsd  (%rax),%xmm0
   14ef8:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   14efc:	48 8d 50 08          	lea    0x8(%rax),%rdx
   14f00:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
   14f04:	f2 0f 5e 45 f0       	divsd  -0x10(%rbp),%xmm0
   14f09:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   14f0d:	83 45 e4 01          	addl   $0x1,-0x1c(%rbp)
   14f11:	8b 45 e4             	mov    -0x1c(%rbp),%eax
   14f14:	3b 45 d4             	cmp    -0x2c(%rbp),%eax
   14f17:	7c d7                	jl     14ef0 <_Z7l2_normPKdiPPd+0x62>
   14f19:	90                   	nop
   14f1a:	90                   	nop
   14f1b:	c9                   	leave  
   14f1c:	c3                   	ret    

0000000000014f1d <_Z7vec_addPKdS0_iPPd>:
   14f1d:	f3 0f 1e fa          	endbr64 
   14f21:	55                   	push   %rbp
   14f22:	48 89 e5             	mov    %rsp,%rbp
   14f25:	48 83 ec 40          	sub    $0x40,%rsp
   14f29:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   14f2d:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
   14f31:	89 55 cc             	mov    %edx,-0x34(%rbp)
   14f34:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
   14f38:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   14f3c:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   14f40:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   14f44:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   14f48:	8b 45 cc             	mov    -0x34(%rbp),%eax
   14f4b:	48 98                	cltq   
   14f4d:	48 c1 e0 03          	shl    $0x3,%rax
   14f51:	48 89 c7             	mov    %rax,%rdi
   14f54:	e8 c7 56 ff ff       	call   a620 <malloc@plt>
   14f59:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   14f5d:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   14f61:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   14f65:	48 89 10             	mov    %rdx,(%rax)
   14f68:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
   14f6f:	eb 38                	jmp    14fa9 <_Z7vec_addPKdS0_iPPd+0x8c>
   14f71:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   14f75:	48 8d 50 08          	lea    0x8(%rax),%rdx
   14f79:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
   14f7d:	f2 0f 10 08          	movsd  (%rax),%xmm1
   14f81:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   14f85:	48 8d 50 08          	lea    0x8(%rax),%rdx
   14f89:	48 89 55 f0          	mov    %rdx,-0x10(%rbp)
   14f8d:	f2 0f 10 00          	movsd  (%rax),%xmm0
   14f91:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   14f95:	48 8d 50 08          	lea    0x8(%rax),%rdx
   14f99:	48 89 55 f8          	mov    %rdx,-0x8(%rbp)
   14f9d:	f2 0f 58 c1          	addsd  %xmm1,%xmm0
   14fa1:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   14fa5:	83 45 e4 01          	addl   $0x1,-0x1c(%rbp)
   14fa9:	8b 45 e4             	mov    -0x1c(%rbp),%eax
   14fac:	3b 45 cc             	cmp    -0x34(%rbp),%eax
   14faf:	7c c0                	jl     14f71 <_Z7vec_addPKdS0_iPPd+0x54>
   14fb1:	90                   	nop
   14fb2:	90                   	nop
   14fb3:	c9                   	leave  
   14fb4:	c3                   	ret    

0000000000014fb5 <_Z7vec_subPKdS0_iPPd>:
   14fb5:	f3 0f 1e fa          	endbr64 
   14fb9:	55                   	push   %rbp
   14fba:	48 89 e5             	mov    %rsp,%rbp
   14fbd:	48 83 ec 40          	sub    $0x40,%rsp
   14fc1:	48 89 7d d8          	mov    %rdi,-0x28(%rbp)
   14fc5:	48 89 75 d0          	mov    %rsi,-0x30(%rbp)
   14fc9:	89 55 cc             	mov    %edx,-0x34(%rbp)
   14fcc:	48 89 4d c0          	mov    %rcx,-0x40(%rbp)
   14fd0:	48 8b 45 d8          	mov    -0x28(%rbp),%rax
   14fd4:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   14fd8:	48 8b 45 d0          	mov    -0x30(%rbp),%rax
   14fdc:	48 89 45 f0          	mov    %rax,-0x10(%rbp)
   14fe0:	8b 45 cc             	mov    -0x34(%rbp),%eax
   14fe3:	48 98                	cltq   
   14fe5:	48 c1 e0 03          	shl    $0x3,%rax
   14fe9:	48 89 c7             	mov    %rax,%rdi
   14fec:	e8 2f 56 ff ff       	call   a620 <malloc@plt>
   14ff1:	48 89 45 f8          	mov    %rax,-0x8(%rbp)
   14ff5:	48 8b 45 c0          	mov    -0x40(%rbp),%rax
   14ff9:	48 8b 55 f8          	mov    -0x8(%rbp),%rdx
   14ffd:	48 89 10             	mov    %rdx,(%rax)
   15000:	c7 45 e4 00 00 00 00 	movl   $0x0,-0x1c(%rbp)
   15007:	eb 38                	jmp    15041 <_Z7vec_subPKdS0_iPPd+0x8c>
   15009:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   1500d:	48 8d 50 08          	lea    0x8(%rax),%rdx
   15011:	48 89 55 e8          	mov    %rdx,-0x18(%rbp)
   15015:	f2 0f 10 00          	movsd  (%rax),%xmm0
   15019:	48 8b 45 f0          	mov    -0x10(%rbp),%rax
   1501d:	48 8d 50 08          	lea    0x8(%rax),%rdx
   15021:	48 89 55 f0          	mov    %rdx,-0x10(%rbp)
   15025:	f2 0f 10 08          	movsd  (%rax),%xmm1
   15029:	48 8b 45 f8          	mov    -0x8(%rbp),%rax
   1502d:	48 8d 50 08          	lea    0x8(%rax),%rdx
   15031:	48 89 55 f8          	mov    %rdx,-0x8(%rbp)
   15035:	f2 0f 5c c1          	subsd  %xmm1,%xmm0
   15039:	f2 0f 11 00          	movsd  %xmm0,(%rax)
   1503d:	83 45 e4 01          	addl   $0x1,-0x1c(%rbp)
   15041:	8b 45 e4             	mov    -0x1c(%rbp),%eax
   15044:	3b 45 cc             	cmp    -0x34(%rbp),%eax
   15047:	7c c0                	jl     15009 <_Z7vec_subPKdS0_iPPd+0x54>
   15049:	90                   	nop
   1504a:	90                   	nop
   1504b:	c9                   	leave  
   1504c:	c3                   	ret    

Disassembly of section .fini:

0000000000015050 <_fini>:
   15050:	f3 0f 1e fa          	endbr64 
   15054:	48 83 ec 08          	sub    $0x8,%rsp
   15058:	48 83 c4 08          	add    $0x8,%rsp
   1505c:	c3                   	ret    
